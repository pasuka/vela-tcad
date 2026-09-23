"""Compare fixed-path experiments, preserving failures and referenced QF precision."""
import argparse
import hashlib
import json
import math
from electrothermal_state import read_bound_record
from pathlib import Path

from evidence_paths import candidate_path


def read(path):
    return read_bound_record(path)


def state_delta(a,b):
    x=a['referenced_state_interleaved'];y=b['referenced_state_interleaved']
    if len(x)!=len(y) or len(x)%4:raise ValueError('State layout differs')
    for state in (a,b):
        for key in ('referenced_state_interleaved','electron_qf_reference_V','hole_qf_reference_V'):
            if any(not math.isfinite(v) for v in state[key]):raise ValueError('Nonfinite comparison state')
    if a['potential_origin_V']!=b['potential_origin_V']:raise ValueError('State origins differ')
    maximum=[0.]*4
    for i in range(len(x)//4):
        for k in (0,3):maximum[k]=max(maximum[k],abs(x[4*i+k]-y[4*i+k]))
        for k,key in ((1,'electron_qf_reference_V'),(2,'hole_qf_reference_V')):
            delta=math.fsum((a[key][i],-b[key][i],x[4*i+k],-y[4*i+k]))
            maximum[k]=max(maximum[k],abs(delta))
    return maximum


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--matrix',type=Path,required=True)
    parser.add_argument('--candidate-path-map',nargs=2)
    parser.add_argument('--bias-digits',type=int,choices=(15,),help='Explicit audited voltage serialization; never interpolate states')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    summary=read(args.matrix/'summary.json');comparisons=[]
    sources={str(args.matrix/'summary.json'):hashlib.sha256((args.matrix/'summary.json').read_bytes()).hexdigest()}
    baseline={r['gate_V']:r for r in summary['runs'] if r['variant']=='baseline'}
    def bias_key(value):
        return float(format(value,'.15g')) if args.bias_digits else value
    state_keys=('state_interleaved','referenced_state_interleaved','electron_qf_reference_V','hole_qf_reference_V','diagnostic_predictor_candidates','diagnostic_tangent_predictor')
    for row in summary['runs']:
        if row['variant']=='baseline':continue
        gate=row['gate_V'];base=baseline[gate]
        base_ledger=read(args.matrix/('baseline_vg%d'%gate)/'results/ledger.json')
        ledger=read(args.matrix/('%s_vg%d'%(row['variant'],gate))/'results/ledger.json')
        accepted=[r for r in base_ledger['runs'] if r['gate']['pass_gate']]
        by_bias={bias_key(r['bias_V']):r for r in accepted}
        if len(by_bias)!=len(accepted):raise ValueError('Duplicate serialized reference bias')
        if row['status']==base['status']=='complete':
            expected={bias_key(v) for v in row['targets_V']}
            actual={bias_key(r['bias_V']) for r in ledger['runs'] if r['gate']['pass_gate']}
            if not expected.issubset(set(by_bias)&actual):raise ValueError('Completed runs are missing requested comparison targets')
        differences=[]
        for run in ledger['runs']:
            key=bias_key(run['bias_V'])
            if not run['gate']['pass_gate'] or key not in by_bias:continue
            directories=[candidate_path(v['directory'],args.candidate_path_map) for v in (run,by_bias[key])]
            states=[read(d/'output.json') for d in directories]
            configs=[read(d/'input.json') for d in directories]
            numerical_keys=('diagnostic_density_update_iterations','diagnostic_density_requires_primary_prediction')
            numerical_differences={k:[c.get(k) for c in configs] for k in numerical_keys if configs[0].get(k)!=configs[1].get(k)}
            for cfg in configs:
                for state_key in state_keys:cfg.pop(state_key,None)
                for numerical_key in numerical_keys:cfg.pop(numerical_key,None)
                if args.bias_digits:
                    for boundary in cfg['boundaries']:
                        if boundary['kind']=='neutral_contact':boundary['value']=bias_key(boundary['value'])
            if configs[0]!=configs[1]:raise ValueError('Non-state point configuration differs')
            delta=state_delta(*states)
            for directory in directories:
                path=directory/'output.json';sources[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
            differences.append(dict(bias_V=run['bias_V'],baseline_bias_V=by_bias[bias_key(run['bias_V'])]['bias_V'],max_physical_state_delta_V_V_V_K=delta,exact_state=not any(delta),numerical_control_differences=numerical_differences))
        item=dict(gate_V=gate,variant=row['variant'],status=row['status'],baseline_status=base['status'],
            common_accepted_states=len(differences),requested_targets=len(row['targets_V']),
            same_requested_targets=row['targets_V']==base['targets_V'],complete_comparison=row['status']==base['status']=='complete',
            failed_attempts=row['failed_attempts'],baseline_failed_attempts=base['failed_attempts'],state_comparison=differences,
            all_compared_states_exact=all(d['exact_state'] for d in differences) if differences else False,
            physical_point_configs_equal=True,numerical_controls_equal=not any(d['numerical_control_differences'] for d in differences),voltage_serialization_digits=args.bias_digits,screening_seconds=row['screening_seconds'])
        for metric in ('newton_updates','line_search_trials','assembly_calls','factorizations','wall_seconds','child_cpu_seconds'):
            if metric in row and metric in base:
                item[metric]=dict(baseline=base[metric],candidate=row[metric],ratio=row[metric]/base[metric] if base[metric] else None)
        comparisons.append(item)
    result=dict(scope='Diagnostic common-state comparison only; screening_seconds overlaps assembly_seconds; no new full-curve acceptance',
        comparisons=comparisons,sources_sha256=sources)
    with args.output.open('x',encoding='utf-8') as output:json.dump(result,output,indent=2)
    for row in comparisons:
        print(json.dumps({k:v for k,v in row.items() if k!='state_comparison'}))


if __name__=='__main__':
    main()
