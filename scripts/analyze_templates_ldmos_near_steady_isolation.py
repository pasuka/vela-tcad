"""Audit frozen near-state controls and four-term QF cancellation exposure."""
import argparse
import json
import math
from pathlib import Path

from analyze_templates_ldmos_predictor_study import state_delta
from run_templates_ldmos_density_projection import read, sha
from run_templates_ldmos_electrothermal_curve import state_gate
from run_templates_ldmos_poisson_initialization import total_cost


def cancellation_audit(cfg,state):
    result={}
    for k,key in ((1,'electron_qf_reference_V'),(2,'hole_qf_reference_V')):
        refs=state[key];x=state['referenced_state_interleaved'];count=0;lost=0;worst=None;worst_relative=None
        for edge in cfg['edge_geometry']:
            if edge['transport_weight']==0:continue
            a,b=edge['nodes'];terms=(refs[b],-refs[a],x[4*b+k],-x[4*a+k])
            exact=math.fsum(terms);ordinary=((refs[b]-refs[a])+x[4*b+k])-x[4*a+k]
            error=abs(ordinary-exact)
            if error:
                count+=1
                if ordinary==0. and exact!=0.:lost+=1
                relative=error/abs(exact) if exact else None
                if relative is not None and (worst_relative is None or relative>worst_relative['relative_error']):
                    worst_relative=dict(nodes=[a,b],ordinary_V=ordinary,compensated_V=exact,
                                        absolute_error_V=error,relative_error=relative)
                if worst is None or error>worst['absolute_error_V']:
                    worst=dict(nodes=[a,b],ordinary_V=ordinary,compensated_V=exact,
                               absolute_error_V=error,relative_error=relative,terms=terms)
        result[key]=dict(affected_transport_edges=count,lost_nonzero_differences=lost,
                        worst_absolute=worst,worst_relative=worst_relative)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--matrix',type=Path,required=True);p.add_argument('--local',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    findings=dict(schema='vela.near_steady_findings.v1',matrices=[],runs=[],checks=[],arithmetic=[])
    loaded={}
    for root,followup in ((args.matrix,False),(args.local,True)):
        s=read(root/'summary.json')
        if s['status']!='completed_diagnostic_matrix' or s.get('reference_followup',False)!=followup:
            raise ValueError('Incomplete/wrong matrix')
        source=Path(s['source_directory']);prior=read(source/'summary.json')
        if sha(source/'summary.json')!=s['source_summary_sha256']:raise ValueError('Prior evidence changed')
        if sha(root/'electrothermal_probe.exe')!=s['binary_sha256'] or s['binary_sha256']!=prior['binary_sha256']:
            raise ValueError('Executable differs')
        if s['runtime_sha256']!=prior['runtime_sha256']:raise ValueError('Runtime differs')
        for name,digest in s['runtime_sha256'].items():
            if sha(root/name)!=digest:raise ValueError('Runtime modified')
        variants=('audit_local','r7_local','density_local') if followup else (
            'audit_terminal','audit_rebased','audit_reference','r7','r7_rebased','density','mass_calibration','mass_qf','mass_density')
        if len(s['runs'])!=2*len(variants) or {(r['vg'],r['variant']) for r in s['runs']}!={(g,v) for g in (4,8) for v in variants}:
            raise ValueError('Incomplete controls')
        findings['matrices'].append(dict(path=str(root),summary_sha256=sha(root/'summary.json')))
        for row in s['runs']:
            source_row=next(a for a in s['sources'] if a['vg']==row['vg'])
            for key in ('reference','original_input'):
                if sha(Path(source_row[key+'_path']))!=source_row[key+'_sha256']:raise ValueError('Source file modified')
            directory=root/row['directory']
            if sha(directory/'input.json')!=row['input_sha256'] or sha(directory/'output.json')!=row['output_sha256']:
                raise ValueError('Run changed')
            cfg,data=read(directory/'input.json'),read(directory/'output.json')
            reference=read(Path(source_row['reference_path']))
            gate=state_gate(data,5.333333333333333);delta=state_delta(data,reference)
            qualified=row['exit_code']==0 and gate['pass_gate'] and all(d<=t for d,t in zip(delta,(1e-8,1e-8,1e-8,1e-7)))
            if qualified!=row['qualified'] or gate!=row['gate'] or delta!=row['max_state_delta']:
                raise ValueError('Qualification changed')
            costs=dict(newton_updates=data['newton_updates'],attempts=len(data['history']),
                trials=sum(h['line_search_trials'] for h in data['history']),assemblies=data['performance']['assembly_calls'],
                factorizations=data['performance']['factorizations'],wall_seconds=row['wall_seconds'])
            if any(row[k]!=v for k,v in costs.items()):raise ValueError('Cost changed')
            pt=data.get('pseudo_transient',{})
            if row['variant'] in ('mass_qf','mass_density') and abs(pt['initial_tau_s']/source_row['terminal_tau_s']-1)>1e-12:
                raise ValueError('Mass time not isolated')
            loaded[row['vg'],row['variant']]=(cfg,data,row,source_row)
            findings['runs'].append(dict(vg=row['vg'],variant=row['variant'],qualified=qualified,cost=costs,
                stop=data['diagnostic_stop'],row_max=data['carrier_row_gate']['max_ratio'],blocks=data['electrical_block_gates'],
                state_change=row['state_change'],reference_delta=delta,
                fallbacks=sum(h.get('density_update_fallback',False) for h in data['history']),
                extra_direction_audits=len(data.get('pseudo_direction_audit',{}).get('directions',[])),
                extra_direction_audit_seconds=data.get('pseudo_direction_audit',{}).get('seconds',0.),
                initial_tau_s=pt.get('initial_tau_s'),final_tau_s=pt.get('next_tau_s')))
    for vg in (4,8):
        cfg,terminal,_,sr=loaded[vg,'audit_terminal']
        prior_source=Path(read(args.matrix/'summary.json')['source_directory'])
        prior_row=next(r for r in read(prior_source/'summary.json')['runs'] if r['vg']==vg and r['variant']=='ptc_defect_model_s1')
        original=read(prior_source/prior_row['directory']/'output.json')
        keys=('referenced_state_interleaved','electron_qf_reference_V','hole_qf_reference_V','residual','carrier_row_gate')
        same=all(terminal[k]==original[k] for k in keys)
        if not same:raise ValueError('Zero-update reconstruction differs from frozen terminal')
        r7=loaded[vg,'r7'][1];density=loaded[vg,'density'][1]
        local_audit=loaded[vg,'audit_local'];local=loaded[vg,'r7_local'];mapped=loaded[vg,'density_local']
        # Numerical model/settings may only vary in the documented diagnostics and state representation.
        diagnostic_keys={'diagnostic_newton_max_iterations','diagnostic_pseudo_transient','diagnostic_pseudo_time_scale',
                         'diagnostic_pseudo_acceptance','diagnostic_pseudo_direction_audit','diagnostic_density_projection'}
        state_keys={'state_interleaved','referenced_state_interleaved','electron_qf_reference_V','hole_qf_reference_V'}
        baseline={k:v for k,v in cfg.items() if k not in diagnostic_keys|state_keys}
        for (g,name),(c,d,_,_) in loaded.items():
            if g==vg and {k:v for k,v in c.items() if k not in diagnostic_keys|state_keys}!=baseline:
                raise ValueError('Unexpected physical configuration change')
        findings['checks'].append(dict(vg=vg,terminal_reassembly_exact=same,
            density_and_qf_final_exact=all(r7[k]==density[k] for k in keys),
            local_density_and_qf_final_exact=all(local[1][k]==mapped[1][k] for k in keys),
            local_input_state_change=state_delta(cfg,local_audit[0]),
            local_references_zeroed=[sum(a!=b and b==0. for a,b in zip(cfg[k],local_audit[0][k]))
                for k in ('electron_qf_reference_V','hole_qf_reference_V')],
            charged_prior_plus_local_audit_and_qf=total_cost(sr['prior_charged_cost'],local_audit[2],local[2])))
        for name in ('audit_terminal','audit_local'):
            c,d,_,_=loaded[vg,name]
            findings['arithmetic'].append(dict(vg=vg,variant=name,qf_edge_difference=cancellation_audit(c,d)))
    findings['scope']='Two near-steady states only. Restarts reset scaling/history; mass first tau equals the prior terminal tau. Rebase changes representation, never copies the reference solution. No final-gate relaxation or full-curve acceleration qualification.'
    with args.output.open('x',encoding='utf-8') as stream:json.dump(findings,stream,indent=2,allow_nan=False)
    print(json.dumps(findings['checks']))


if __name__=='__main__':main()
