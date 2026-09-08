"""Independent eight-point qualification, restart checks and field localization."""
import argparse
import ast
import copy
from decimal import Decimal
import math
from pathlib import Path
import re
import subprocess
import numpy as np
import validate_simplemos_joint_eight_point_20260907 as launch

b=launch.b; v=launch.v; a=launch.a; d=launch.d; p=launch.p
LOCAL=launch.LOCAL; OUT=launch.OUT
REPORT=p.REPO/'docs/validation/simplemos_joint_geometry_eight_point_validation_2026-09-07.md'


def cases():
    answer=[]
    for out,root in ((b.OUT,b.LOCAL),(OUT,LOCAL)):
        for original in a.read(out/'contract.json')['cases']:
            c=copy.deepcopy(original); c['root']=str(root/c['key']); c['evidence_root']=str(out)
            c['jobs']=[j for j in c['jobs'] if j['arm'] in ('baseline','joint')]; answer.append(c)
    assert len(answer)==8
    return answer


def support(c):
    geo=d.matrix.spatial.m73.Geometry(c['device'])
    mask=d.matrix.spatial.old.m78.supports(c['device'],geo,.05)[0]['all_si'].copy(); mask[geo.contact_nodes]=False
    assert np.count_nonzero(mask)==907
    return geo,mask


def physical(row,field):
    car={'phin':'electron','phip':'hole'}.get(field)
    if car and car+'_qf_reference_V' in row:
        return Decimal(row[car+'_qf_reference_V'])+Decimal(row[car+'_qf_increment_V'])
    return Decimal(row[field])


def delta_states(first,second,mask):
    diffs={field:max(abs(float(physical(x,field)-physical(y,field))) for x,y,keep in zip(first,second,mask) if keep) for field in ('psi','phin','phip')}
    density=max(abs(float(x[k])/float(y[k])-1) for x,y,keep in zip(first,second,mask) if keep for k in ('electrons_m3','holes_m3'))
    return dict(psi_max_V=diffs['psi'],phin_max_V=diffs['phin'],phip_max_V=diffs['phip'],density_max_relative=density)


def stripped(cfg):
    cfg=copy.deepcopy(cfg)
    for k in ('state_file','output_state_file'):cfg.pop(k)
    cfg['solver']['local_update_diagnostics'].pop('csv_file')
    return cfg


def qualify(c,dest):
    geo,mask=support(c); s=a.read(dest/'config.status.json'); q=a.read(dest/'all_row.status.json')
    assert stripped(a.read(dest/'config.json'))==stripped(a.read(Path(c['base'])/'config.json'))
    post=a.read(dest/'all_row.json'); assert post['solver']['global_continuity_closure']==dict(mode='enforce',tolerance=1e-6,source_floor=1e-10)
    rows=[]
    for i,(r,keep) in enumerate(zip(d.ordered(dest/'all_row.csv',geo.count),mask)):
        if not keep:continue
        for car in ('electron','hole'):
            scale=max(float(r[car+'_flux_abs_sum']),abs(float(r[car+'_recombination'])),abs(float(r[car+'_impact'])))
            ratio=abs(float(r[car+'_residual']))/scale if scale else math.inf
            rows.append(dict(node_id=i,carrier=car,x_um=float(r['x']),y_um=float(r['y']),row_ratio=ratio,
                residual=float(r[car+'_residual']),absolute_edge_flow=float(r[car+'_flux_abs_sum']),net_edge_flow=float(r[car+'_flux']),
                SRH=float(r[car+'_recombination']),impact=float(r[car+'_impact']),scale=scale))
    assert len(rows)==q['carrier_row_convergence']['qualified_row_count']==1814
    bad=sum(r['row_ratio']>1e-6 for r in rows); zero=sum(r['scale']==0 for r in rows)
    assert bad==q['carrier_row_convergence']['violation_count']
    cc=s['contact_currents_A_per_um']; current=cc['drain']; kcl=abs(math.fsum(cc.values()))/abs(current)
    ok=s['exit_code']==q['exit_code']==0 and s['converged'] and bad==zero==0 and q['global_continuity_closure']['satisfied'] and kcl<=1e-8
    info=dict(current_A_per_um=current,qualified=ok,active_rows=1814,row_violations=bad,zero_scale_rows=zero,
              max_row_ratio=max(r['row_ratio'] for r in rows),kcl_over_Id=kcl,iterations=s['iterations'])
    return info,sorted(rows,key=lambda r:r['row_ratio'],reverse=True)


def restart():
    for f in (OUT/'freeze.json',OUT/'preflight_evidence.json',b.OUT/'validation_evidence.json'):a.verify(f)
    jobs=[]; files=[Path(__file__).resolve(),OUT/'freeze.json',b.OUT/'validation_evidence.json']
    for c in cases():
        for label in ('replacement','replacement_native'):
            old=Path(c['root'])/'joint'/label; info,rows=qualify(c,old)
            if info['max_row_ratio']<5e-7:continue
            assert info['qualified']
            dest=LOCAL/'near_threshold_restart'/c['key']/label
            cfg=a.read(old/'config.json'); cfg.update(state_file=str(old/'state.csv'),output_state_file=str(dest/'state.csv'))
            cfg['solver']['local_update_diagnostics']['csv_file']=str(dest/'updates.csv'); a.write(dest/'config.json',cfg)
            post=a.read(old/'all_row.json'); post.update(state_file=str(dest/'state.csv'),output_csv=str(dest/'all_row.csv')); a.write(dest/'all_row.json',post)
            jobs.append(dict(case=c,label=label,old=str(old),dest=str(dest),original_info=info,worst_original_row=rows[0]))
            files += [old/n for n in ('state.csv','config.status.json','all_row.csv','all_row.status.json')]+[dest/'config.json',dest/'all_row.json']
    a.write(OUT/'restart_contract.json',dict(status='frozen_before_restart',trigger_ratio=5e-7,jobs=jobs,unchanged_gates=True))
    files.append(OUT/'restart_contract.json'); d.matrix.freeze(OUT/'restart_freeze.json',files)
    result=[]
    for j in jobs:
        c=j['case']; dest=Path(j['dest']); old=Path(j['old']); v.execute(dest/'config.json',c,1.,1.); v.execute(dest/'all_row.json',c,1.,1.)
        info,ranked=qualify(c,dest); geo,mask=support(c)
        changes=delta_states(d.ordered(dest/'state.csv',geo.count),d.ordered(old/'state.csv',geo.count),mask)
        iddiff=abs(info['current_A_per_um']/j['original_info']['current_A_per_um']-1)
        ok=info['qualified'] and max(changes[k] for k in ('psi_max_V','phin_max_V','phip_max_V'))<=1e-6 and changes['density_max_relative']<=1e-4 and iddiff<=1e-6
        result.append(dict(key=c['key'],label=j['label'],original_max_row_ratio=j['original_info']['max_row_ratio'],
            restarted_max_row_ratio=info['max_row_ratio'],restarted_worst_node=ranked[0]['node_id'],restarted_worst_carrier=ranked[0]['carrier'],
            Id_relative=iddiff,**changes,restart_iterations=info['iterations'],qualified=ok))
        print('Unchanged-gate restart',c['key'],j['label'],'qualified',ok,'max row',info['max_row_ratio'],flush=True)
    a.write_csv(OUT/'near_threshold_restarts.csv',result)
    files += [x for x in (LOCAL/'near_threshold_restart').rglob('*') if x.is_file()]+[OUT/'near_threshold_restarts.csv']
    d.matrix.freeze(OUT/'restart_evidence.json',files)


def analyze():
    for f in (OUT/'freeze.json',OUT/'preflight_evidence.json',OUT/'restart_evidence.json',b.OUT/'validation_evidence.json'):a.verify(f)
    for device in ('n19','n23'):
        for name in ('edges.txt','nodes.txt'):assert a.sha(b.old.LOCAL/'ratios'/device/name)==a.sha(b.LOCAL/'ratios'/device/name)
    checks=[]; rankedrows=[]; comparisons=[]; fd=[]; invariance=[]; fields=[]; peaks=[]; jvp=[]; local=[]
    for c in cases():
        geo,mask=support(c); root=Path(c['root']); mapped=d.ordered(Path(c['native_initial']),geo.count)
        weights=geo.volumes['signed_si'][mask]; weights=weights/np.sum(weights); values={}; passed={}
        for job in c['jobs']:
            dest=Path(job['config']).parent; info,ranked=qualify(c,dest)
            meta=dict(key=c['key'],device=c['device'],vg=c['vg'],vd=c['vd'],arm=job['arm'],label=job['label'])
            checks.append(dict(**meta,**info)); rankedrows += [dict(**meta,rank=i+1,**r) for i,r in enumerate(ranked[:10])]
            values[(job['arm'],job['label'])]=info['current_A_per_um']; passed[(job['arm'],job['label'])]=info['qualified']
            saved=next(x for x in a.rows(Path(c['evidence_root'])/'dc.csv') if x['key']==c['key'] and x['arm']==job['arm'] and x['label']==job['label'])
            assert info['qualified']==(saved['qualified']=='True') and info['current_A_per_um']==float(saved['current_A_per_um'])
            linear=a.rows(dest/'linear_summary.csv')
            for it in {r['iteration'] for r in linear}:assert [int(r['correction']) for r in linear if r['iteration']==it]==list(range(5))
            if job['label'] not in ('zero','replacement','replacement_native'):continue
            fun=a.read(dest/'functional.status.json'); assert fun['exit_code']==0
            porterr=abs(fun['current_A_per_um']/fun['contact_current_extractor_A_per_um']-1)
            jp=launch.check.jvp_metrics(dest/'jvp.csv',c['key'],job['arm']+'/'+job['label'])[0]; jvp+=jp
            comparisons.append(dict(**meta,initialization=job['initialization'],**info,native_Id_A_per_um=c['native_Id_A_per_um'],
                signed_Id_relative_error=info['current_A_per_um']/c['native_Id_A_per_um']-1,port_functional_relative=porterr,
                jvp_changed_blocks_qualified=all(r['qualified'] for r in jp)))
            assert porterr<=1e-8
            if job['initialization']!='vela':continue
            state=d.ordered(dest/'state.csv',geo.count); terms=d.ordered(dest/'all_row.csv',geo.count)
            for field in ('psi','phin','phip'):
                delta=np.array([float(physical(x,field)-physical(y,field)) for x,y in zip(state,mapped)])
                order=sorted(np.flatnonzero(mask),key=lambda i:abs(delta[i]),reverse=True)
                fields.append(dict(**meta,field=field,qualified_state=info['qualified'],free_Si_max_absolute_V=abs(delta[order[0]]),
                    max_node_id=int(order[0]),Si_volume_weighted_RMS_V=float(np.sqrt(np.sum(weights*delta[mask]**2)))))
                for i in list(dict.fromkeys(order[:5]+[320,324,338])):
                    r=terms[i]; s=state[i]; n=mapped[i]
                    record=dict(**meta,field=field,node_id=int(i),x_um=float(r['x']),y_um=float(r['y']),delta_V=delta[i],
                        Vela_physical_V=float(physical(s,field)),native_physical_V=float(physical(n,field)),
                        Vela_electrons_m3=float(s['electrons_m3']),Vela_holes_m3=float(s['holes_m3']),
                        native_coherent_electrons_m3=float(n['electrons_m3']),native_coherent_holes_m3=float(n['holes_m3']))
                    for car in ('electron','hole'):
                        scale=max(float(r[car+'_flux_abs_sum']),abs(float(r[car+'_recombination'])),abs(float(r[car+'_impact'])))
                        record.update({car+'_row_ratio':abs(float(r[car+'_residual']))/scale,car+'_absolute_edge_flow':float(r[car+'_flux_abs_sum']),
                            car+'_net_edge_flow':float(r[car+'_flux']),car+'_SRH':float(r[car+'_recombination'])})
                    if i in order[:5]:peaks.append(dict(rank=order.index(i)+1,**record))
                    if i in (320,324,338):local.append(record)
        center=values[('baseline','zero')]; drift=abs(center-c['base_Id_A_per_um']); assert abs(center/c['base_Id_A_per_um']-1)<=1e-8
        pred=float(next(x['unit_prediction_A_per_um'] for x in a.rows(Path(c['evidence_root'])/'predictions.csv') if x['key']==c['key'] and x['arm']=='joint'))
        odd={n:(values[('joint','plus_'+n)]-values[('joint','minus_'+n)])/2 for n in ('full','half')}; linearity=abs(odd['full']/(2*odd['half'])-1)
        for name,alpha in [('full',.001),('half',.0005)]:
            plus,minus=values[('joint','plus_'+name)],values[('joint','minus_'+name)]
            err=abs(odd[name]/(alpha*pred)-1); even=abs(((plus+minus)/2-center)/odd[name]); snr=abs(odd[name])/max(drift,1e-300); dd=abs(math.log10(center/c['base_Id_A_per_um']))
            signs=(plus-center)*pred>0 and (minus-center)*pred<0
            ok=passed[('baseline','zero')] and all(passed[('joint',s+'_'+name)] for s in ('plus','minus')) and err<=.001 and linearity<=.001 and even<=.01 and snr>=100 and dd<=1e-5 and signs
            fd.append(dict(key=c['key'],vg=c['vg'],amplitude=name,prediction_relative_error=err,two_amplitude_relative=linearity,even_over_odd=even,
                signal_over_zero_drift=snr,zero_drift_dex=dd,signs_correct=signs,unit_prediction_A_per_um=pred,actual_odd_A_per_um=odd[name],qualified=ok))
        change=delta_states(d.ordered(root/'joint/replacement/state.csv',geo.count),d.ordered(root/'joint/replacement_native/state.csv',geo.count),mask)
        iddiff=abs(values[('joint','replacement')]/values[('joint','replacement_native')]-1)
        ok=passed[('joint','replacement')] and passed[('joint','replacement_native')] and max(change[k] for k in ('psi_max_V','phin_max_V','phip_max_V'))<=1e-6 and change['density_max_relative']<=1e-4 and iddiff<=1e-6
        invariance.append(dict(key=c['key'],vg=c['vg'],**change,Id_relative=iddiff,qualified=ok))
        print('Independently checked seven states, response, initialization and fields',c['key'],flush=True)
    pairs=[]
    for vg in (.8,1.):
        for vd in (.05,1.):
            vals={arm:{dev:next(x for x in comparisons if x['vg']==vg and x['vd']==vd and x['device']==dev and x['arm']==arm and x['initialization']=='vela') for dev in ('n19','n23')} for arm in ('baseline','joint')}
            errs={arm:math.log10(val['n23']['current_A_per_um']/val['n19']['current_A_per_um'])-math.log10(val['n23']['native_Id_A_per_um']/val['n19']['native_Id_A_per_um']) for arm,val in vals.items()}
            hi=abs(vals['joint']['n23']['signed_Id_relative_error'])<abs(vals['baseline']['n23']['signed_Id_relative_error'])
            lo=abs(vals['joint']['n19']['signed_Id_relative_error'])<=abs(vals['baseline']['n19']['signed_Id_relative_error'])
            better=abs(errs['joint'])<abs(errs['baseline'])
            pairs.append(dict(vg=vg,vd=vd,baseline_pair_error_dex=errs['baseline'],joint_pair_error_dex=errs['joint'],high_NWell_improved=hi,low_NWell_not_worse=lo,pair_improved=better,qualified=hi and lo and better))
    restarts=a.rows(OUT/'near_threshold_restarts.csv')
    qualified=all(x['qualified'] for group in (checks,fd,invariance,pairs,jvp) for x in group) and all(x['qualified']=='True' for x in restarts)
    for name,rows in [('independent_qualification',checks),('ranked_carrier_rows',rankedrows),('eight_point_comparison',comparisons),('response_calibration',fd),('initialization_invariance',invariance),('field_comparison',fields),('field_peaks',peaks),('local_fields',local),('replacement_jvp',jvp),('NWell_pairing',pairs)]:a.write_csv(OUT/(name+'.csv'),rows)
    metrics=dict(new_DC=28,new_qualified_DC=sum(x['qualified'] for x in checks if x['vg']==.8),combined_joint_study_DC=len(checks),combined_qualified_DC=sum(x['qualified'] for x in checks),
        combined_independently_checked_carrier_rows=len(checks)*1814,total_amplitudes=len(fd),qualified_amplitudes=sum(x['qualified'] for x in fd),
        initialization_pairs=len(invariance),qualified_initialization_pairs=sum(x['qualified'] for x in invariance),restart_count=len(restarts),qualified_restarts=sum(x['qualified']=='True' for x in restarts),
        max_row_ratio=max(x['max_row_ratio'] for x in checks),max_kcl_over_Id=max(x['kcl_over_Id'] for x in checks),
        max_response_prediction_error=max(x['prediction_relative_error'] for x in fd),max_response_two_amplitude_difference=max(x['two_amplitude_relative'] for x in fd),
        max_initialization_phi_V=max(max(x[k] for k in ('psi_max_V','phin_max_V','phip_max_V')) for x in invariance),
        max_initialization_density_relative=max(x['density_max_relative'] for x in invariance),max_initialization_Id_relative=max(x['Id_relative'] for x in invariance),
        eight_point_joint_gate_passed=qualified,M82_released=False,M83_released=False)
    a.write(OUT/'analysis_evidence.json',dict(status='completed_bounded_eight_point_extension',metrics=metrics,
        input_hashes={a.rel(f):a.sha(f) for f in [Path(__file__).resolve()]+[x for root in (LOCAL,OUT) for x in root.rglob('*') if x.is_file()]}))
    print(metrics,flush=True)


def seal():
    for file in ('freeze.json','preflight_evidence.json','restart_evidence.json','analysis_evidence.json'):a.verify(OUT/file)
    a.verify(b.OUT/'validation_evidence.json')
    for script in (Path(__file__).resolve(),Path(launch.__file__).resolve()):ast.parse(script.read_text(encoding='utf-8'))
    for target in re.findall(r'\]\(([^)]+)\)',REPORT.read_text(encoding='utf-8')):
        if '://' not in target:
            path=(REPORT.parent/target.split('#')[0]).resolve()
            if path.name!='validation_evidence.json':assert path.is_file(),path
    diff=subprocess.run(['git','-c','core.fsmonitor=false','diff','--numstat','--ignore-space-at-eol'],cwd=p.REPO,capture_output=True,text=True)
    assert diff.returncode==0 and diff.stdout.strip()=='127\t1\tsrc/tools/vela_example_runner.cpp'
    files=[REPORT,Path(__file__).resolve(),Path(launch.__file__).resolve(),b.OUT/'validation_evidence.json']
    files += [x for root in (LOCAL,OUT) for x in root.rglob('*') if x.is_file()]
    a.write(OUT/'validation_evidence.json',dict(status='completed_bounded_eight_point_joint_geometry_validation',date='2026-09-07',
        metrics=a.read(OUT/'analysis_evidence.json')['metrics'],production_changes=False,acceptance_changes=False,
        input_hashes={a.rel(f):a.sha(f) for f in sorted(set(files))},
        limitations=['Vg=0.8/1 V only; not a complete 0-1 V curve.','No new native simulation; reuse sealed qualified native exports.',
            'Stable isolated analytical SG psi derivative reused; no production changes.','No full native cell mobility, dielectric K, or SRH source volume replacement.',
            'Directional JVP checks are not exhaustive arbitrary-column Jacobian qualification.','M82/M83 not released.']))
    a.verify(OUT/'validation_evidence.json'); print('Eight-point evidence sealed, old evidence unchanged.',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('action',choices=('restart','analyze','seal')); globals()[parser.parse_args().action]()
