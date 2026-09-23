"""Summarize completed paired state-I/O controls without pooling executables."""
import argparse
import json
from pathlib import Path


def analyze(root):
    batch=json.loads((root/'summary.json').read_text(encoding='utf-8'))
    if batch['status']!='completed' or len(batch['cases'])!=batch['planned_curves']:
        raise ValueError('Incomplete paired batch')
    if not batch['paired_timing']:raise ValueError('Not a paired timing experiment')
    cases=batch['cases'];control=batch['control'];candidate=batch['candidate'];pairs=[]
    keys={(c['repeat'],c['gate']) for c in cases}
    for repeat,gate in sorted(keys):
        selected={c['mode']:c for c in cases if (c['repeat'],c['gate'])==(repeat,gate)}
        a,b=selected[control],selected[candidate]
        for c in (a,b):
            if c['status']!='completed' or not c['audit']['integrity_pass'] or not c['trajectory_exact']:
                raise ValueError('State/trajectory audit failed')
            if c['audit']['exact_points']!=batch['points'] or c['max_potential_difference_V']>1e-8:
                raise ValueError('Incomplete or changed state')
        ca,cb=a['performance']['counters'],b['performance']['counters']
        differences={k:[ca.get(k),cb.get(k)] for k in sorted(set(ca)|set(cb)) if ca.get(k)!=cb.get(k)}
        pairs.append(dict(repeat=repeat,gate_V=gate,control_seconds=a['wall_seconds'],candidate_seconds=b['wall_seconds'],
            saving_percent=100*(1-b['wall_seconds']/a['wall_seconds']),
            control_updates=a['audit']['Newton_updates'],candidate_updates=b['audit']['Newton_updates'],
            services=[a['audit']['services'],b['audit']['services']],counter_differences=differences,
            state_difference_V=[a['max_potential_difference_V'],b['max_potential_difference_V']],
            worker_roundtrip_seconds=[c['parent_timing']['worker_roundtrip']['inclusive'] for c in (a,b)],
            parent_nonworker_seconds=[c['wall_seconds']-c['parent_timing']['worker_roundtrip']['inclusive'] for c in (a,b)],
            memory_seed_stats=b.get('memory_seed_stats')))
    a=sum(p['control_seconds'] for p in pairs);b=sum(p['candidate_seconds'] for p in pairs)
    return dict(control=control,candidate=candidate,points_per_curve=batch['points'],curves=len(cases),
        repeats=len({c['repeat'] for c in cases}),single_round_is_not_stability=True,
        combined_saving_percent=100*(1-b/a),runtime_sha256=batch['runtime_sha256'],pairs=pairs)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--batch',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();result=analyze(a.batch)
    a.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='runtime_sha256'},indent=2))
