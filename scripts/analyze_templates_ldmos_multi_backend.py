"""Analyze arbitrary backend names without conflating solver and ordering suffixes."""
import argparse
import statistics
from pathlib import Path
from run_templates_ldmos_linked_d5 import read,write
from analyze_templates_ldmos_isothermal_backends import repeat_statistics,factor_statistics

def matrix_summary(root):
    summary=read(root/'summary.json');groups={}
    for case in summary['cases']:
        name=case['backend'];group=groups.setdefault(name,dict(samples=[],failures=[]))
        if case['status']!='pass':group['failures'].append(case);continue
        rows=read(root/f'r{case["round"]}_{name}.json')['systems'];times={}
        for r in rows:
            for stage in r['profiling']['stages']:
                times[stage['name']]=times.get(stage['name'],0.)+stage['total_ns']/1e9
        peaks=[r['profiling'].get('resources',{}).get('process_peak_working_set_bytes') for r in rows]
        peaks=[p for p in peaks if p is not None]
        group['samples'].append(dict(round=case['round'],wall_seconds=case['wall_seconds'],seconds=times,
            process_cpu_seconds=case.get('cpu_seconds'),
            matrices=len(rows),analyses=rows[-1]['cumulative_analyses'],
            worst_solver_input_backward_error=max(r['quality']['scaled']['normwise_backward_error'] for r in rows),
            raw_available=all(r['quality'].get('raw_available',True) for r in rows),
            peak_process_bytes=max(peaks) if peaks else None))
    for group in groups.values():
        samples=group['samples']
        if not samples:continue
        times=[r['seconds']['linear.total'] for r in samples]
        group['median_seconds']={key:statistics.median(r['seconds'][key] for r in samples) for key in samples[0]['seconds']}
        group['service_range_over_mean_percent']=100*(max(times)-min(times))/statistics.mean(times)
    return dict(status=summary['status'],threads=summary['threads'],groups=groups)

def curve_summary(root):
    summary=read(root/'summary.json');cases=[]
    for case in summary['cases']:
        if case.get('status')!='pass' or 'audit' not in case:continue
        child=root/case['name'];ledger=read(child/'fixed/ledger.json')
        profiles=[read(child/r['case']/'performance_profile.json') for r in ledger['runs']]
        stages={}
        for profile in profiles:
            for stage in profile['stages']:stages[stage['name']]=stages.get(stage['name'],0.)+stage['total_ns']/1e9
        prefix=case['name'].split('_',3)
        cases.append(dict(name=case['name'],backend=case['backend'],profile=prefix[0],gate=int(prefix[1][2:]),round=int(prefix[2][1:]),
            solver_threads=case.get('solver_threads'),blas_threads=case.get('blas_threads'),
            wall_seconds=case['wall_seconds'],parent_cpu_seconds=case['cpu_seconds']['total'],
            child_cpu_seconds=case['audit']['child_cpu_seconds'],
            process_tree_cpu_seconds=case['cpu_seconds']['total']+case['audit']['child_cpu_seconds'],audit=case['audit'],
            Newton_updates=ledger['total_Newton_updates'],rollbacks=len(ledger['rollbacks']),point_services=len(ledger['runs']),
            stages=stages,factor_statistics=factor_statistics(profiles),system_cpu_samples=case.get('system_cpu_samples',[])))
    pairs=[]
    for candidate in cases:
        controls=[r for r in cases if r['backend']=='umfpack' and all(r[k]==candidate[k] for k in ['profile','gate','round'])]
        if controls and candidate['backend']!='umfpack':
            control=controls[0]
            pairs.append(dict(profile=candidate['profile'],gate=candidate['gate'],round=candidate['round'],candidate=candidate['backend'],
                wall_saving_vs_umfpack_percent=100*(1-candidate['wall_seconds']/control['wall_seconds']),
                cpu_saving_vs_umfpack_percent=100*(1-candidate['process_tree_cpu_seconds']/control['process_tree_cpu_seconds']),
                newton_change=candidate['Newton_updates']-control['Newton_updates']))
    return dict(status=summary['status'],cases=cases,pairs=pairs,repeats=repeat_statistics(root,cases),
                comparisons=summary['comparisons'],joint_qualifications=summary.get('joint_qualifications',{}))

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path);p.add_argument('--matrix',action='store_true');a=p.parse_args()
    result=matrix_summary(a.root) if a.matrix else curve_summary(a.root)
    write(a.root/'multi_backend_analysis.json',result)
    print('ANALYSIS_COMPLETE',a.root)

if __name__=='__main__':main()
