"""Check full-trajectory switching, frozen controls, qualification and total costs."""
import argparse
import copy
import json
from pathlib import Path
from run_templates_ldmos_density_projection import read,sha
from run_templates_ldmos_poisson_initialization import total_cost
from run_templates_ldmos_electrothermal_curve import state_gate
from analyze_templates_ldmos_predictor_study import state_delta


STATE_RESULT=('state_interleaved','referenced_state_interleaved','electron_qf_reference_V',
              'hole_qf_reference_V','residual','carrier_row_gate','electrical_block_gates','diagnostic_stop','newton_updates')


def strip_switch(history):
    result=copy.deepcopy(history)
    for h in result:h.pop('near_steady_qf_active',None)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--matrix',type=Path,required=True)
    p.add_argument('--prepared-control',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    root=args.matrix;s=read(root/'summary.json')
    if s['status']!='completed_diagnostic_matrix':raise ValueError('Incomplete matrix')
    if sha(root/'electrothermal_probe.exe')!=s['binary_sha256']:raise ValueError('Executable changed')
    for name,h in s['runtime_sha256'].items():
        if sha(root/name)!=h:raise ValueError('Runtime changed')
    oldroot=Path(s['control_directory']);old=read(oldroot/'summary.json')
    prepared=read(args.prepared_control/'summary.json')
    if sha(oldroot/'summary.json')!=s['control_summary_sha256']:raise ValueError('Control changed')
    if old['status']!='completed_diagnostic_matrix' or prepared['status']!='completed_diagnostic_matrix':
        raise ValueError('Incomplete previous controls')
    if len(s['runs'])!=6 or {(r['vg'],r['variant']) for r in s['runs']}!={(g,v) for g in (4,8) for v in ('baseline','switch_raw','switch_prepared')}:
        raise ValueError('Missing candidate runs')
    def checked(base,row,name=None):
        directory=base/(name or row['directory'])
        if sha(directory/'input.json')!=row['input_sha256'] or sha(directory/'output.json')!=row['output_sha256']:
            raise ValueError('Input/output changed')
        return read(directory/'input.json'),read(directory/'output.json')
    findings=dict(schema='vela.near_switch_findings.v1',summary_sha256=sha(root/'summary.json'),
        prepared_control_summary_sha256=sha(args.prepared_control/'summary.json'),runs=[])
    for row in s['runs']:
        cfg,data=checked(root,row);vg=row['vg'];name='vg%d_p%d'%(vg,11 if vg==4 else 9)
        if sha(Path(row['reference_path']))!=row['reference_sha256']:raise ValueError('Reference changed')
        gate=state_gate(data,5.333333333333333);delta=state_delta(data,read(Path(row['reference_path'])))
        qualified=row['exit_code']==0 and gate['pass_gate'] and all(d<=t for d,t in zip(delta,(1e-8,1e-8,1e-8,1e-7)))
        if qualified!=row['qualified'] or delta!=row['max_state_delta'] or gate!=row['gate']:raise ValueError('Qualification differs')
        counts=dict(newton_updates=data['newton_updates'],attempts=len(data['history']),
            trials=sum(h['line_search_trials'] for h in data['history']),assemblies=data['performance']['assembly_calls'],
            factorizations=data['performance']['factorizations'],wall_seconds=row['wall_seconds'])
        if any(row[k]!=v for k,v in counts.items()):raise ValueError('Cost differs')
        previous_variant='baseline' if row['variant']=='baseline' else 'ptc_defect_model_s1'
        previous_root=args.prepared_control if row['variant']=='switch_prepared' else oldroot
        previous_summary=prepared if row['variant']=='switch_prepared' else old
        previous_row=next(r for r in previous_summary['runs'] if r['vg']==vg and r['variant']==previous_variant)
        pcfg,prior=checked(previous_root,previous_row,previous_variant+'_'+name)
        compare=dict(cfg);compare.pop('diagnostic_near_steady_qf_switch',None)
        if compare!=pcfg:raise ValueError('Unexpected configuration difference')
        switch=data.get('near_steady_qf_switch',{});events=switch.get('events',[])
        if len(events)>1:raise ValueError('Repeated switch')
        prefix=data['history'] if not events else [h for h in data['history'] if h['iteration']<events[0]['next_iteration']]
        if strip_switch(prefix)!=prior['history'][:len(prefix)]:raise ValueError('Pre-switch history changed')
        if not events:
            if any(data[k]!=prior[k] for k in STATE_RESULT):raise ValueError('Inactive switch changed trajectory')
        else:
            event=events[0]
            if event['merit_before']>=1e-9:raise ValueError('Early switch')
            if not all(b['satisfied'] for b in prior['history'][event['next_iteration']-1]['iteration_trace']['before_gates']['blocks']):
                raise ValueError('Switch before original block gates passed')
            for h in data['history']:
                active=h['iteration']>=event['next_iteration']
                if h['near_steady_qf_active']!=active:raise ValueError('Invalid mode history')
                if active and ('pseudo_transient' in h or 'projection_trials' in h or h.get('density_update_attempted',False)):
                    raise ValueError('Mass or density survived switch')
            if any(v['iteration']>=event['next_iteration'] for v in data['pseudo_transient']['steps']):
                raise ValueError('Mass continued after switch')
        prep_cost=[]
        if row['variant']=='switch_prepared':
            prep=next(a for a in s['preparations'] if a['vg']==vg)
            for key in ('preparation','audit'):checked(root,prep[key]);prep_cost.append(prep[key])
            if not prep['gate']['pass_gate']:raise ValueError('Preparation failed')
        charged=total_cost(*prep_cost,row)
        if charged!=row['charged_cost']:raise ValueError('Preparation cost omitted')
        findings['runs'].append(dict(vg=vg,variant=row['variant'],qualified=qualified,cost=counts,charged_cost=charged,
            switch=switch,pre_switch_history_exact=True,untouched_trajectory_exact=not bool(events),
            row_max=data['carrier_row_gate']['max_ratio'],gate=gate,max_state_delta=delta,
            stop=data['diagnostic_stop'],post_switch_updates=sum(h['accepted'] for h in data['history'] if h.get('near_steady_qf_active',False)),
            steady_ratio=(data['pseudo_transient']['steps'][-1]['fixed_residual_after']/data['pseudo_transient']['steps'][0]['fixed_residual_before']
                          if data.get('pseudo_transient',{}).get('steps') else None)))
    findings['scope']='Complete original-seed attempts, with Poisson preparation charged where used. Optional solver transition only; original steady gates. Diagnostic single-run timing is not a repeated full-curve benchmark.'
    with args.output.open('x',encoding='utf-8') as stream:json.dump(findings,stream,indent=2,allow_nan=False)
    print(json.dumps([dict(vg=r['vg'],variant=r['variant'],qualified=r['qualified'],cost=r['charged_cost'],switch=r['switch']) for r in findings['runs']]))


if __name__=='__main__':main()
