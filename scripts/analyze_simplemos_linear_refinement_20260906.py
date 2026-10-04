"""Keep linear accuracy, DC qualification, current changes and seed invariance distinct."""
from pathlib import Path
from decimal import Decimal
import itertools
import math
import numpy as np
import validate_simplemos_linear_refinement_20260906 as m

a=m.a;d=m.d


def field_difference(left,right,mask):
    result={}
    result['psi_max_V']=max(abs(float(l['psi'])-float(r['psi'])) for l,r,keep in zip(left,right,mask) if keep)
    for carrier in ('electron','hole'):
        def qf(row):return Decimal(row[carrier+'_qf_reference_V'])+Decimal(row[carrier+'_qf_increment_V'])
        result[carrier+'_qf_max_V']=float(max(abs(qf(l)-qf(r)) for l,r,keep in zip(left,right,mask) if keep))
        key='electrons_m3' if carrier=='electron' else 'holes_m3'
        result[carrier+'_density_max_symmetric_relative']=max(abs(float(l[key])-float(r[key]))/max(abs(float(l[key])),abs(float(r[key])),1e-300) for l,r,keep in zip(left,right,mask) if keep)
    return result


def main():
    a.verify(m.OUT/'freeze.json');a.verify(m.v.OUT/'validation_evidence.json')
    rows=[];linear=[];trials=[];groups={};files=[Path(__file__).resolve(),m.OUT/'freeze.json',m.OUT/'runs.csv',m.OUT/'disabled_sentinels.csv']
    for job in a.read(m.OUT/'contract.json')['jobs']:
        path=Path(job['config']).parent;source=Path(job['source']);s=a.read(path/'config.status.json');r=a.read(path/'result.json')
        old=a.read(source/'result.json');os=a.read(source/'config.status.json');gate=a.read(path/'all_row.status.json')['carrier_row_convergence']
        geo=d.matrix.spatial.m73.Geometry(job['device']);mask=d.matrix.spatial.old.m78.supports(job['device'],geo,.05)[0]['all_si'].copy();mask[geo.contact_nodes]=False
        state=d.ordered(path/'state.csv',geo.count);before=d.ordered(source/'state.csv',geo.count);terms=d.ordered(path/'all_row.csv',geo.count)
        zero=sum(max(abs(float(t[c+'_flux_abs_sum'])),abs(float(t[c+'_recombination'])),abs(float(t[c+'_impact'])))==0 for t,keep in zip(terms,mask) if keep for c in ('electron','hole'))
        assert gate['qualified_row_count']==2*int(sum(mask))==1814
        delta=(r['current_A_per_um']-old['current_A_per_um'])/max(abs(old['current_A_per_um']),1e-300)
        row=dict(tag=job['tag'],device=job['device'],vd=job['vd'],vg=job['vg'],initialization=job['initialization'],
            baseline_converged=os['converged'],refined_converged=s['converged'],baseline_qualified=old['all_row_qualified'],refined_qualified=r['all_row_qualified'] and zero==0,
            baseline_iterations=os['iterations'],refined_iterations=s['iterations'],baseline_failure=os['failure_reason'],refined_failure=s['failure_reason'],
            baseline_violations=old['all_row_violations'],refined_violations=r['all_row_violations'],baseline_max_ratio=old['all_row_max_ratio'],refined_max_ratio=r['all_row_max_ratio'],
            electron_violations=sum(v['carrier']=='electron' for v in gate['violations']),hole_violations=sum(v['carrier']=='hole' for v in gate['violations']),
            zero_scale_rows=zero,baseline_Id_A_per_um=old['current_A_per_um'],refined_Id_A_per_um=r['current_A_per_um'],Id_change_over_baseline_abs=delta,
            baseline_kcl_over_Id=old['kcl_over_Id'],refined_kcl_over_Id=r['kcl_over_Id'],all_row_global_satisfied=r['all_row_global_satisfied'],
            **field_difference(before,state,mask))
        history=s['newton_trace'];summary=a.rows(path/'linear_summary.csv');iters=sorted({int(x['iteration']) for x in summary})
        final_linear=[]
        for it in iters:
            g=[x for x in summary if int(x['iteration'])==it];assert [int(x['correction']) for x in g]==list(range(5))
            assert all(int(x['active_rows'])==1814 for x in g)
            first,last=g[0],g[-1];t=next(t for t in history if t['iteration']==it)
            entry=dict(tag=job['tag'],iteration=it,initial_linear_max=float(first['max_carrier_linear_ratio']),final_linear_max=float(last['max_carrier_linear_ratio']),
                initial_backward_max=float(first['componentwise_backward_max']),final_backward_max=float(last['componentwise_backward_max']),
                final_linear_rows_above_1e_minus_8=int(last['rows_above_1e_minus_8']),zero_scale_rows=int(last['zero_scale_rows']),
                linear_monitor_passed=int(last['zero_scale_rows'])==0 and float(last['max_carrier_linear_ratio'])<=1e-8,
                step_inf_V=float(last['step_inf_V']),accepted=t['accepted'],event=t['event'])
            final_linear.append(entry);linear.append(entry)
        row.update(linear_steps=len(final_linear),linear_monitor_passed_steps=sum(x['linear_monitor_passed'] for x in final_linear),
            max_final_linear_ratio=max(x['final_linear_max'] for x in final_linear),max_final_backward=max(x['final_backward_max'] for x in final_linear))
        last=max(history,key=lambda t:t['iteration']);prev=max((t for t in history if t['iteration']<last['iteration']),key=lambda t:t['iteration'])
        scales=np.array([max(history[0]['blocks'][key],1.) for key in ('psi','phin','phip')]);blocks=np.array([prev['blocks'][key] for key in ('psi','phin','phip')])
        merit=float(np.linalg.norm(blocks/scales));row['last_poisson_merit_fraction']=float((blocks[0]/scales[0])**2)/max(merit**2,1e-300)
        for t in last['trials']:trials.append(dict(tag=job['tag'],iteration=last['iteration'],baseline_merit=merit,**t))
        rows.append(row);groups.setdefault((job['case'],job['vg']),[]).append((row,state,mask))
        files += [path/n for n in ('config.status.json','state.csv','all_row.csv','all_row.status.json','legacy.status.json','result.json','linear_summary.csv','linear_selected.csv','updates.csv')]
        print(job['tag'],old['all_row_violations'],'->',row['refined_violations'],'qualified',row['refined_qualified'],'linear',row['linear_monitor_passed_steps'],'/',row['linear_steps'],flush=True)
    pairs=[]
    for (case,vg),group in groups.items():
        full=len(group)==4 and {x[0]['initialization'] for x in group}=={'vela','native','hole_plus','hole_minus'}
        qualified=full and all(x[0]['refined_qualified'] for x in group)
        for left,right in itertools.combinations(group,2):
            l,sl,mask=left;r,sr,_=right;fields=field_difference(sl,sr,mask)
            current=abs(l['refined_Id_A_per_um']-r['refined_Id_A_per_um'])/max(abs(l['refined_Id_A_per_um']),abs(r['refined_Id_A_per_um']),1e-300)
            numeric=max(fields['psi_max_V'],fields['electron_qf_max_V'],fields['hole_qf_max_V'])<=1e-6 and max(fields['electron_density_max_symmetric_relative'],fields['hole_density_max_symmetric_relative'])<=1e-4 and current<=1e-6
            pairs.append(dict(case=case,vg=vg,left=l['initialization'],right=r['initialization'],full_four_seeds=full,all_four_qualified=qualified,
                both_states_qualified=l['refined_qualified'] and r['refined_qualified'],Id_symmetric_relative=current,**fields,numerical_invariance=numeric,invariance_qualified=qualified and numeric))
    for name,data in (('comparison',rows),('linear_iterations',linear),('last_trials',trials),('initialization_pairs',pairs)):a.write_csv(m.OUT/(name+'.csv'),data)
    a.write(m.OUT/'result.json',dict(input_hashes={a.rel(p):a.sha(p) for p in sorted(set(files))},completed_states=len(rows),
        baseline_qualified=sum(r['baseline_qualified'] for r in rows),refined_qualified=sum(r['refined_qualified'] for r in rows),
        newly_qualified=sum(r['refined_qualified'] and not r['baseline_qualified'] for r in rows),lost_qualification=sum(r['baseline_qualified'] and not r['refined_qualified'] for r in rows),
        linear_steps=len(linear),linear_monitor_passed_steps=sum(r['linear_monitor_passed'] for r in linear),
        qualified_four_seed_cases=sum(len(g)==4 and all(r[0]['refined_qualified'] for r in g) and all(p['numerical_invariance'] for p in pairs if p['case']==case and p['vg']==vg) for (case,vg),g in groups.items()),
        production_changes=False,old_gate_changed=False,remote_runs=0,
        limitation='Only seven selected initial states at three bias points; failed state current agreement does not qualify a physical comparison or resolve Id-Vg discrepancy.'))


if __name__=='__main__':main()
