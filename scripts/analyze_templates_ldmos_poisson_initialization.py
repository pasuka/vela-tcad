"""Verify preparation isolation, original gates, PTC decisions and charged costs."""
import argparse
import json
import math
from pathlib import Path

from analyze_templates_ldmos_pseudo_transient import checked_matrix, read, sha
from analyze_templates_ldmos_predictor_study import state_delta
from run_templates_ldmos_electrothermal_curve import state_gate
from run_templates_ldmos_poisson_initialization import (
    VARIANTS, preparation_gate, restored_config, total_cost, variant_config)


def metrics(row, data):
    steps=data.get('pseudo_transient',{}).get('steps',[])
    trials=[t for s in steps for t in s.get('defect_trials',[])]
    for s in steps:
        for t in s.get('defect_trials',[]):
            if not all(math.isfinite(t[k]) for k in ('alpha','defect_norm','model_error','steady_fixed_norm')):
                raise ValueError('Nonfinite trial metric')
            if t['accepted'] != (t['defect_norm'] <= (1.-1e-4*t['alpha'])*s['fixed_residual_before']):
                raise ValueError('Recorded defect acceptance differs from Armijo')
    cost=dict(newton_updates=data['newton_updates'],attempts=len(data['history']),
              trials=sum(h['line_search_trials'] for h in data['history']),
              assemblies=data['performance']['assembly_calls'],
              factorizations=data['performance']['factorizations'],wall_seconds=row['wall_seconds'])
    closed=[h['iteration'] for h in data['history']
            if h.get('iteration_trace',{}).get('after_gates',{}).get('blocks')
            and all(b['satisfied'] for b in h['iteration_trace']['after_gates']['blocks'])]
    return dict(cost=cost,stop=data['diagnostic_stop'],
        final_blocks=data['electrical_block_gates'],final_row_gate=data['carrier_row_gate'],
        final_residual_blocks=data['residual_blocks'],
        accepted_defect_steps=sum(s.get('defect_used',False) for s in steps),
        accepted_steady_increases=sum(s.get('defect_used',False) and s['fixed_residual_after']>s['fixed_residual_before'] for s in steps),
        defect_trial_checks=len(trials),fallbacks=sum(h.get('density_update_fallback',False) for h in data['history']),
        floor_attempts=sum(h['scaled_l2_before']<1e-9 for h in data['history']),
        first_electrical_blocks_pass_iteration=closed[0] if closed else None,
        fixed_residual_final_ratio=steps[-1]['fixed_residual_after']/steps[0]['fixed_residual_before'] if steps else None,
        initial_tau_s=data.get('pseudo_transient',{}).get('initial_tau_s'),
        final_tau_s=data.get('pseudo_transient',{}).get('next_tau_s'),
        first_step=steps[0] if steps else None)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--matrix',type=Path,required=True)
    parser.add_argument('--manifest',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();root=args.matrix
    summary=read(root/'summary.json')
    if summary['schema']!='vela.poisson_initialization.v1' or summary['status']!='completed_diagnostic_matrix':
        raise ValueError('Incomplete or wrong experiment')
    control_root=Path(summary['control_directory']);old,controls=checked_matrix(control_root)
    if sha(control_root/'summary.json')!=summary['control_summary_sha256']:
        raise ValueError('Control summary changed')
    if sha(args.manifest)!=summary['manifest_sha256'] or old['manifest_sha256']!=summary['manifest_sha256']:
        raise ValueError('Manifest changed')
    if sha(root/'electrothermal_probe.exe')!=summary['binary_sha256'] or summary['binary_sha256']!=old['binary_sha256']:
        raise ValueError('Executable differs')
    if summary['runtime_sha256']!=old['runtime_sha256']:
        raise ValueError('Runtime differs')
    for name,digest in summary['runtime_sha256'].items():
        if sha(root/name)!=digest:raise ValueError('Runtime changed')
    cases={c['gate_V']:c for c in read(args.manifest) if c['variant']=='linear_baseline'
           and (c['gate_V'],c['baseline_index']) in ((4,11),(8,9))}
    if len(summary['preparations'])!=2 or {p['vg'] for p in summary['preparations']}!={4,8}:
        raise ValueError('Missing preparation')
    def checked(row):
        directory=root/row['directory']
        if sha(directory/'input.json')!=row['input_sha256'] or sha(directory/'output.json')!=row['output_sha256']:
            raise ValueError('Recorded point changed')
        cfg,data=read(directory/'input.json'),read(directory/'output.json')
        m=metrics(row,data)
        if any(row[k]!=v for k,v in m['cost'].items()):raise ValueError('Recorded cost differs')
        return cfg,data,m
    def instrument(cfg):return dict(cfg,performance_profiling=True,diagnostic_iteration_trace=True)
    result=dict(schema='vela.poisson_initialization_findings.v1',summary_sha256=sha(root/'summary.json'),preparations=[],comparisons=[])
    prepared_by_gate={}
    for item in summary['preparations']:
        vg=item['vg'];original=read(Path(cases[vg]['input']))
        if sha(Path(cases[vg]['input']))!=item['source_input_sha256']:raise ValueError('Source changed')
        cfg,prepared,pm=checked(item['preparation'])
        if cfg!=instrument(dict(original,solve_mode='poisson',initialization='provided_state')):
            raise ValueError('Unexpected preparation configuration change')
        acfg,audit,am=checked(item['audit'])
        restored=restored_config(original,prepared)
        if acfg!=instrument(dict(restored,diagnostic_newton_max_iterations=0)):
            raise ValueError('Original coupled model not restored')
        gate=preparation_gate(original,prepared,audit)
        if gate!=item['gate'] or item['ready']!=gate['pass_gate']:
            raise ValueError('Preparation gate differs')
        prepared_by_gate[vg]=(restored,item)
        result['preparations'].append(dict(vg=vg,gate=gate,cost=total_cost(item['preparation'],item['audit']),
            coupled_blocks=audit['electrical_block_gates'],coupled_row_gate=audit['carrier_row_gate'],
            coupled_residual_blocks=audit['residual_blocks'],coupled_merit=audit.get('merit')))
    expected={(v,p['vg']) for p in summary['preparations'] if p['ready'] for v in VARIANTS}
    if len(summary['runs'])!=len(expected) or {(r['variant'],r['vg']) for r in summary['runs']}!=expected:
        raise ValueError('Incomplete candidate matrix')
    control_map={(r['variant'],r['vg']):(r,d) for r,d in controls}
    for row in summary['runs']:
        cfg,data,new=checked(row);vg=row['vg'];restored,item=prepared_by_gate[vg]
        if cfg!=instrument(variant_config(restored,row['variant'])):
            raise ValueError('Unexpected candidate configuration change')
        gate=state_gate(data,cases[vg]['bias_V']);delta=state_delta(data,read(Path(cases[vg]['expected'])))
        qualified=row['exit_code']==0 and gate['pass_gate'] and all(d<=t for d,t in zip(delta,(1e-8,1e-8,1e-8,1e-7)))
        if row['qualified']!=qualified or row['gate']!=gate or row['max_state_delta']!=delta:
            raise ValueError('Final qualification differs')
        charged=total_cost(item['preparation'],item['audit'],row)
        if charged!=row['charged_cost']:raise ValueError('Preparation cost omitted')
        cr,cd=control_map[row['variant'],vg]
        result['comparisons'].append(dict(variant=row['variant'],vg=vg,qualified=qualified,
            control_qualified=cr['qualified'],control=metrics(cr,cd),prepared=new,charged_cost=charged,
            final_gate=gate,max_state_delta=delta))
    result['scope']='Poisson-only consistency at fixed QF/T, not complete algebraic consistency; fixed residual scales reset per solve. Every variant is charged full preparation and audit cost.'
    with args.output.open('x',encoding='utf-8') as stream:json.dump(result,stream,indent=2,allow_nan=False)
    print(json.dumps([dict(variant=c['variant'],vg=c['vg'],qualified=c['qualified'],
         updates=c['prepared']['cost']['newton_updates'],charged=c['charged_cost']['newton_updates']) for c in result['comparisons']]))


if __name__=='__main__':main()
