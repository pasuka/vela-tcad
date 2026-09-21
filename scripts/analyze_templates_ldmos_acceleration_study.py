"""Summarize all passes/failures; timings exclude capture reading and error checks."""
import argparse
import json
import math
from pathlib import Path
import statistics
import struct


def summarize(root):
    report=json.loads((root/'summary.json').read_text())
    rows=[]
    for config in report['configurations']:
        cases=[c for c in report['cases'] if c['configuration']['name']==config['name']]
        details=[]
        for c in cases:
            path=root/(c['name']+'.json')
            if path.is_file(): details.append(json.loads(path.read_text()))
        passed=[d for d in details if d['status']=='pass']
        row=dict(configuration=config['name'],runs=len(cases),passed=len(passed),
                 returncodes=[c['returncode'] for c in cases])
        if passed:
            row['linear_seconds_median']=statistics.median(d['linear_seconds'] for d in passed)
            row['linear_seconds_range']=[min(d['linear_seconds'] for d in passed),max(d['linear_seconds'] for d in passed)]
            row['factor_seconds_median']=statistics.median(sum(s['factor_seconds'] for s in d['systems']) for d in passed)
            row['factor_entries_sum_median']=statistics.median(sum(s['factor_entries'] for s in d['systems']) for d in passed)
            for key in ('preparation_seconds','analysis_seconds','reorder_seconds','reorder_other_seconds'):
                if all(key in s for d in passed for s in d['systems']):
                    row[key+'_median']=statistics.median(sum(s[key] for s in d['systems']) for d in passed)
            if all('vendor_phase_seconds' in s for d in passed for s in d['systems']):
                phases=set().union(*(s['vendor_phase_seconds'].keys() for d in passed for s in d['systems']))
                row['vendor_phase_seconds_median']={k:statistics.median(sum(s['vendor_phase_seconds'].get(k,0.) for s in d['systems']) for d in passed) for k in sorted(phases)}
            row['analysis_counts']=[sum(not s['analysis_reused'] for s in d['systems']) for d in passed]
            if all('factor_fill_ratio' in s for d in passed for s in d['systems']):
                row['factor_fill_ratio_median']=statistics.median(s['factor_fill_ratio'] for d in passed for s in d['systems'])
            row['krylov_iterations_two_rhs_sum']=[sum(s['two_rhs_krylov_iterations'] for s in d['systems']) for d in passed]
            row['cpu_seconds_median']=statistics.median(c['cpu_seconds']['total'] for c in cases if c['returncode']==0)
            row['process_wall_seconds_median']=statistics.median(c['wall_seconds'] for c in cases if c['returncode']==0)
            row['max_backward_error']=max(s['quality']['scaled']['normwise_backward_error'] for d in passed for s in d['systems'])
            # For solver-input-only captures, report absolute L2 error bounds as
            # well: relative errors of round-off-sized directions are misleading.
            max_bound=0.
            for d in passed:
                for s in d['systems']:
                    if s['quality']['raw_available']: continue
                    data=Path(s['input']).read_bytes(); n=struct.unpack_from('<Q',data,24)[0]
                    reference=struct.unpack_from(f'<{n}d',data,len(data)-8*n)
                    for block,ratio in enumerate(s['quality']['relative_solver_step_difference_by_block']):
                        v=reference[block*(n//3):(block+1)*(n//3)]
                        max_bound=max(max_bound,ratio*math.sqrt(sum(x*x for x in v)))
            row['max_solver_step_block_l2_difference']=max_bound
        rows.append(row)
    return dict(scope=report['scope'],status=report['status'],configurations=rows)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path);a=p.parse_args()
    result=summarize(a.root)
    (a.root/'analysis.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    for row in result['configurations']: print(json.dumps(row))
