"""Streaming row/step diagnostics for the isolated T470p controls."""
import argparse
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path


def read(path):return json.loads(path.read_text(encoding='utf-8'))
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def rows(path):
    with path.open(newline='',encoding='utf-8') as f:yield from csv.DictReader(f)


def row_groups(path):
    for key,group in itertools.groupby(rows(path),key=lambda r:(r['event'],r['iteration'])):
        yield key,list(group)


def diagnostic(case,weights):
    dest=Path(case['directory']);src=Path(case['source'])
    terminal={r['node_id']:r for r in rows(src/'state.csv')}
    cfg=read(dest/'control.json')
    if cfg.get('scaling',{}).get('mode')!='unit_scaling':
        raise ValueError('This frozen diagnostic requires the TCAD concentration conversion')
    # Newton row_trace incorrectly labels internal cm^-3 as m^-3. Unlike the
    # restart writer it does not call internalConcentrationToM3. Keep raw files
    # intact and explicitly convert; audit against this run's terminal state.
    density_factor=1e6
    # Membership is checked against actual row gate below, not inferred from density.
    groups=[];phase=-1
    for (event,it),group in row_groups(dest/'row_trace.csv'):
        if event=='initial':phase+=1
        result=dict(event=event,phase=phase,iteration=int(it),fixed_carrier_l2=0.,worst_rows=[])
        squared=[];worst=[]
        for r in group:
            node=r['node_id']
            for carrier,density in (('electron','n_m3'),('hole','p_m3')):
                residual=float(r[carrier+'_residual']);scale=float(r[carrier+'_scale'])
                srh=float(r[carrier+'_srh']);impact=float(r[carrier+'_impact'])
                fixed=max(weights[(node,carrier)],1e-30)
                squared.append((residual/fixed)**2)
                # This frozen deck qualifies by nonzero source; no density/flux exemptions.
                if srh==0 and impact==0:continue
                ratio=float(r[carrier+'_ratio'])
                worst.append(dict(node=int(node),carrier=carrier,ratio=ratio,residual=residual,scale=scale,
                    density_m3=float(r[density])*density_factor,flux=float(r[carrier+'_flux']),srh=srh,impact=impact,
                    cancellation_ratio=abs(residual)/max(abs(float(r[carrier+'_flux'])),abs(srh),abs(impact),1e-300)))
        result['fixed_carrier_l2']=math.sqrt(math.fsum(squared))
        result['worst_rows']=sorted(worst,key=lambda r:r['ratio'],reverse=True)[:3]
        if event=='initial':
            result['initial_density_error']={}
            for col,density in (('electrons_m3','n_m3'),('holes_m3','p_m3')):
                diffs=[(math.log10(float(r[density])*density_factor)-math.log10(float(terminal[r['node_id']][col])),int(r['node_id']))
                    for r in group if float(r[density])>0 and float(terminal[r['node_id']][col])>0]
                result['initial_density_error'][col]=dict(log10_rms=math.sqrt(math.fsum(d*d for d,n in diffs)/len(diffs)),
                    worst_log10=max(diffs,key=lambda x:abs(x[0])),nodes=len(diffs))
        groups.append(result)
    own_terminal={r['node_id']:r for r in rows(dest/'state.csv')}
    unit_check=max(abs(float(r[d])*density_factor-float(own_terminal[r['node_id']][c]))/
                   max(abs(float(own_terminal[r['node_id']][c])),1.) for r in group
                   for c,d in (('electrons_m3','n_m3'),('holes_m3','p_m3')))
    if unit_check>1e-9:raise ValueError('Terminal density unit audit failed')
    trace=list(rows(dest/'newton_iterations.csv'))
    gate=cfg['solver']['block_absolute_convergence']
    qualified=[r for r in trace if all(float(r[k])<=gate[g] for k,g in (
        ('block_psi','psi_residual_ceiling'),('block_phin','electron_residual_ceiling'),('block_phip','hole_residual_ceiling')))
        and float(r['carrier_row_max_ratio'])>1e-8]
    local=list(rows(dest/'local_updates.csv')) if (dest/'local_updates.csv').exists() else []
    local_groups=[];local_phase=0;last_iteration=0
    for it,rs in itertools.groupby(local,key=lambda r:r['iteration']):
        rs=list(rs)
        if int(it)<last_iteration:local_phase+=1
        last_iteration=int(it)
        local_groups.append(dict(iteration=int(it),phase=local_phase,rows=len(rs),
            max_raw_step_V=max(abs(float(r['raw_linear_step_V'])) for r in rs),
            max_applied_step_V=max(abs(float(r['applied_step_V'])) for r in rs),
            capped_rows=sum(abs(float(r['raw_linear_step_V'])-float(r['capped_step_V']))>1e-12 for r in rs),
            max_componentwise_backward_error=max(float(r['raw_linear_componentwise_backward_error']) for r in rs)))
    curve=list(rows(dest/'curve.csv'))[0]
    profile=read(dest/'performance_profile.json')
    return dict(directory=str(dest),gate=case['gate'],target_V=case['target_V'],variant=case['variant'],
        density_raw_to_m3=density_factor,terminal_density_conversion_relative_error=unit_check,
        internal_recovery={k:v for k,v in curve.items() if k.startswith('carrier_row_recovery')},
        profile_counters=profile['counters'],
        updates=case['updates'],row_tail=case.get('row_only_tail_updates',0),
        trace_sha256=sha(dest/'row_trace.csv'),qualified_block_states=qualified,
        groups=groups,selected_node_linear_updates=local_groups)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--steps',type=Path,required=True);p.add_argument('--diagnostics',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();s=read(a.steps/'summary.json');d=read(a.diagnostics/'summary.json')
    result=dict(steps=[],diagnostics=[],input_hashes={str(a.steps/'summary.json'):sha(a.steps/'summary.json'),
                                                  str(a.diagnostics/'summary.json'):sha(a.diagnostics/'summary.json')})
    for case in s['cases']:
        item={k:v for k,v in case.items() if k!='runs'}
        curve_rows=[r for run in case['runs'] for r in rows(Path(run['directory'])/'curve.csv')]
        item.update(services=len(case['runs']),updates=sum(r['updates'] for r in case['runs']),
            failed_attempts=sum(r['failed_attempts'] for r in case['runs']),
            wall_seconds=sum(r['wall_seconds'] for r in case['runs']),
            row_tail=sum(r.get('row_only_tail_updates',0) for r in case['runs']),
            internal_recovery_points=sum(int(r.get('carrier_row_recovery_attempted','0')) for r in curve_rows),
            internal_density_passes=sum(int(r.get('carrier_row_recovery_density_passes','0')) for r in curve_rows))
        result['steps'].append(item)
    weights={}
    for case in d['cases']:
        key=(case['gate'],case['target_V'])
        if key not in weights:
            if case['variant']!='guarded':raise ValueError('Guarded common scale must be first')
            first=next(row_groups(Path(case['directory'])/'row_trace.csv'))[1]
            weights[key]={(r['node_id'],c):float(r[c+'_scale']) for r in first for c in ('electron','hole')}
        result['diagnostics'].append(diagnostic(case,weights[key]))
        print('analyzed',case['gate'],case['target_V'],case['variant'],flush=True)
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False),encoding='utf-8')


if __name__=='__main__':main()
