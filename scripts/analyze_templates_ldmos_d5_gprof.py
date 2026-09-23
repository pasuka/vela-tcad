"""Cross-check gprof self samples against uninstrumented full-curve stage timers."""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path


def read(path):return json.loads(path.read_text(encoding='utf-8'))
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def aggregate(directory):
    ledger=read(directory/'fixed/ledger.json')
    if ledger['status']!='completed':raise ValueError('Incomplete curve')
    counters=Counter();seconds=Counter();calls=Counter();blas=set();peak=0
    for run in ledger['runs']:
        data=read(directory/run['case']/'performance_profile.json')
        counters.update(data['counters'])
        for stage in data['stages']:
            seconds[stage['name']]+=stage['total_ns']/1e9
            calls[stage['name']]+=stage['calls']
        obs=data['observations'].get('linear.openblas_threads')
        if obs:blas.update((obs['min'],obs['max']))
        peak=max(peak,data['resources'].get('process_peak_working_set_bytes',0))
    if blas!={1}:raise ValueError('BLAS thread observation differs from one')
    return dict(counters=dict(counters),stage_seconds=dict(seconds),stage_calls=dict(calls),
                child_wall_seconds=ledger['total_child_wall_seconds'],
                child_cpu_seconds=ledger['total_child_cpu_seconds'],
                process_peak_working_set_bytes=peak,blas_threads=sorted(blas),
                ledger_sha256=sha(directory/'fixed/ledger.json'))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--profiles',type=Path,required=True)
    p.add_argument('--baseline',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();profile=read(a.profiles/'summary.json');baseline=read(a.baseline/'summary.json')
    if profile['status']!='completed' or baseline['status']!='completed':raise ValueError('Incomplete batch')
    report=dict(status='pass',profiles_summary_sha256=sha(a.profiles/'summary.json'),
                baseline_summary_sha256=sha(a.baseline/'summary.json'),gates={})
    keys=('newton.updates','dd.jacobian_calls','dd.residual_calls','linear.factorize_calls',
          'linear.analyze_calls','linear.solve_calls','newton.line_search_attempts')
    for gate in (4,8):
        case=next(c for c in profile['cases'] if c['gate']==gate)
        if case['status']!='completed':raise ValueError('Incomplete profile case')
        directory=Path(case['directory'])
        instrumented=aggregate(directory);ordinary=aggregate(a.baseline/f'baseline_vg{gate}')
        mismatch={k:[ordinary['counters'].get(k,0),instrumented['counters'].get(k,0)] for k in keys
                  if ordinary['counters'].get(k,0)!=instrumented['counters'].get(k,0)}
        if mismatch:raise ValueError('Instrumented numerical work differs: '+str(mismatch))
        with (directory/'gprof_hotspots.csv').open(newline='',encoding='utf-8') as f:
            samples=[dict(r,percent_time=float(r['percent_time']),self_seconds=float(r['self_seconds'])) for r in csv.DictReader(f)]
        samples.sort(key=lambda r:r['self_seconds'],reverse=True)
        total=sum(r['self_seconds'] for r in samples)
        if total<=0:raise ValueError('No gprof samples')
        oldcase=next(c for c in baseline['cases'] if c['gate']==gate and c['mode']=='baseline')
        row=dict(ordinary=ordinary,instrumented=instrumented,sampled_self_seconds=total,
                 ordinary_wall_seconds=oldcase['wall_seconds'],instrumented_wall_seconds=case['instrumented_wall_seconds'],
                 max_potential_difference_V=case['max_potential_difference_V'],trajectory_exact=case['trajectory_exact'],
                 samples=samples,top_samples=samples[:30],numerical_counters_equal=True,
                 gmon_sha256=sha(directory/'gmon.out'),flat_sha256=sha(directory/'gprof_flat.txt'))
        row['instrumentation_self_seconds']=sum(r['self_seconds'] for r in samples
                                               if any(k in r['function'] for k in ('mcount','__fentry__','__monstartup')))
        report['gates'][str(gate)]=row
    a.output.write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps({k:dict(sampled_self_seconds=v['sampled_self_seconds'],
                            instrumentation_self_seconds=v['instrumentation_self_seconds'],
                            ordinary_wall_seconds=v['ordinary_wall_seconds'],
                            instrumented_wall_seconds=v['instrumented_wall_seconds']) for k,v in report['gates'].items()}),flush=True)


if __name__=='__main__':main()
