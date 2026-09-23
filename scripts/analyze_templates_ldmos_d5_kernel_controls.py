"""Validate completed kernel controls and report same-binary timing pairs."""
import argparse
import hashlib
import json
import math
from pathlib import Path


def analyze(summary):
    if summary['status'] != 'completed':
        raise ValueError('Incomplete batch')
    cases = summary['cases']
    if len(cases) != summary['planned_curves']:
        raise ValueError('Missing planned curves')
    keys = [(c['repeat'], c['mode'], c['gate']) for c in cases]
    if len(set(keys)) != len(keys):
        raise ValueError('Duplicate curve')
    rows = []
    for c in cases:
        p = c['performance']; count = p['counters']; stage = p['stage_seconds']
        if (c['status'] != 'completed' or not c['trajectory_exact'] or
            not math.isfinite(c['max_potential_difference_V']) or c['max_potential_difference_V'] > 1e-8 or
            not c['audit']['integrity_pass'] or c['exact_points'] != summary['points'] or
            p['blas_threads'] != [1]):
            raise ValueError('Unqualified curve')
        if summary['points'] == 31 and not all(v['pass_all'] for v in c['verdicts'].values()):
            raise ValueError('Original full-curve gates failed')
        hits = count.get('jacobian.endpoint_cache_hits', 0)
        misses = count.get('jacobian.endpoint_cache_misses', 0)
        rows.append(dict(repeat=c['repeat'], mode=c['mode'], gate=c['gate'],
                         wall_seconds=c['wall_seconds'],
                         child_wall_seconds=p['child_wall_seconds'], child_cpu_seconds=p['child_cpu_seconds'],
                         updates=c['audit']['Newton_updates'],
                         pattern_builds=count.get('jacobian.pattern_build_calls', 0),
                         structure_hits=count.get('jacobian.structure_cache_hits', 0),
                         pattern_seconds=stage.get('jacobian.pattern_build', 0),
                         structure_check_seconds=stage.get('jacobian.structure_cache_check', 0),
                         endpoint_misses=misses, endpoint_hit_fraction=hits/max(1,hits+misses),
                         edge_physics_seconds=stage.get('jacobian.edge_physics', 0),
                         jacobian_seconds=stage.get('dd.jacobian', stage.get('jacobian', 0)),
                         peak_memory_bytes=p['process_peak_working_set_bytes']))
    pairs=[]
    if summary.get('paired_timing',True):
        control=summary.get('control','baseline');candidate=summary.get('candidate','structure')
        if control==candidate:raise ValueError('Control and candidate must differ')
        for repeat in sorted({c['repeat'] for c in cases}):
            for gate in (4,8):
                a=next(r for r in rows if (r['repeat'],r['mode'],r['gate'])==(repeat,control,gate))
                b=next(r for r in rows if (r['repeat'],r['mode'],r['gate'])==(repeat,candidate,gate))
                if a['updates']!=b['updates']:raise ValueError('Numerical work changed')
                if any(not math.isfinite(t) or t<=0 for t in (a['wall_seconds'],b['wall_seconds'])):
                    raise ValueError('Invalid timing')
                pairs.append(dict(repeat=repeat,gate=gate,control=control,candidate=candidate,
                                  control_wall_seconds=a['wall_seconds'],candidate_wall_seconds=b['wall_seconds'],
                                  wall_reduction_percent=100*(1-b['wall_seconds']/a['wall_seconds'])))
    if summary['points']==31:
        expected=len(cases)//2
        if len(summary['joint'])!=expected or not all(all(j['checks'].values()) for j in summary['joint']):
            raise ValueError('Missing or failed dual-gate joint gates')
    return dict(status='pass',points=summary['points'],rows=rows,pairs=pairs,
                paired_timing=summary.get('paired_timing',True),
                qualification_only=not summary.get('paired_timing',True),
                caveat='One pair per gate is screening evidence, not stable acceleration.')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('summary',type=Path);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();raw=a.summary.read_bytes();result=analyze(json.loads(raw))
    result['summary_sha256']=hashlib.sha256(raw).hexdigest()
    a.output.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
