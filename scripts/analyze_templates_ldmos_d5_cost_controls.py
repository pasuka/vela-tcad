"""Summarize complete cost-continuation controls, retaining failed work."""
import argparse
import hashlib
import json
import math
from pathlib import Path
from ldmos_cost_continuation import work_units


def analyze(batch):
    if batch['status']!='completed' or len(batch['cases'])!=batch['planned_curves']:
        raise ValueError('Incomplete batch')
    rows=[]
    for c in batch['cases']:
        if c['status']!='completed' or not c['audit']['integrity_pass']:
            raise ValueError('Unqualified curve')
        if not math.isfinite(c['max_potential_difference_V']) or c['max_potential_difference_V']>1e-8:
            raise ValueError('State equivalence failed')
        if c['exact_points']!=batch['points']:raise ValueError('Missing exact point')
        if c['mode'] in ('baseline','files') and not c['trajectory_exact']:raise ValueError('Changed baseline')
        if c['mode'] in ('files','efficiency') and not c.get('cold_final_audit'):
            raise ValueError('Missing uncached final audit')
        p=c['performance'];t=c['parent_timing']
        if p['blas_threads']!=[1]:raise ValueError('Unexpected BLAS threads')
        total=t['controller']['inclusive'];exclusive=sum(x['exclusive'] for x in t.values())
        if not all(math.isfinite(x[k]) and x[k]>=-1e-9 for x in t.values() for k in ('inclusive','exclusive')):
            raise ValueError('Invalid timing')
        if abs(total-exclusive)>1e-5:raise ValueError('Timing coverage overlaps or has gaps')
        if batch['points']==31 and not all(v['pass_all'] for v in c['verdicts'].values()):
            raise ValueError('Original curve gates failed')
        rows.append(dict(gate=c['gate'],mode=c['mode'],repeat=c['repeat'],wall_seconds=c['wall_seconds'],
                         updates=c['audit']['Newton_updates'],services=c['audit']['services'],
                         advances=c['audit']['advances'],rollbacks=c['audit']['rollbacks'],
                         work_units=work_units(p['counters']),max_state_V=c['max_potential_difference_V'],
                         factorization_calls=p['counters'].get('linear.factorize_calls',0),
                         jacobian_calls=p['counters'].get('dd.jacobian_calls',0),
                         residual_calls=p['counters'].get('dd.residual_calls',0),
                         child_wall_seconds=p['child_wall_seconds'],child_cpu_seconds=p['child_cpu_seconds'],
                         controller_seconds=total,worker_roundtrip_seconds=t['worker_roundtrip']['inclusive'],
                         parent_exclusive_seconds=total-t['worker_roundtrip']['inclusive'],
                         file_cache_stats=c.get('file_cache_stats',{}),
                         parent_breakdown={k:v for k,v in t.items() if k!='worker_roundtrip'}))
    if batch['points']==31 and (len(batch['joint'])!=len(rows)//2 or not all(all(x['checks'].values()) for x in batch['joint'])):
        raise ValueError('Missing joint acceptance')
    return dict(status='pass',points=batch['points'],rows=rows)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('summary',type=Path);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();raw=a.summary.read_bytes();r=analyze(json.loads(raw));r['summary_sha256']=hashlib.sha256(raw).hexdigest()
    a.output.write_text(json.dumps(r,indent=2),encoding='utf-8');print(json.dumps(r,indent=2))

if __name__=='__main__':main()
