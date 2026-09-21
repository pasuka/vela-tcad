"""Inspect actual permutations and repeated timings without conflating LU policies."""
import argparse
import json
from pathlib import Path
import statistics
from run_templates_ldmos_linked_d5 import read,write

def summarize(root):
    summary=read(root/'summary.json');groups={}
    for c in summary['cases']:
        if c['returncode']!=0:continue
        data=read(root/f"r{c['round']}_{c['mode']}.json");groups.setdefault(c['mode'],[]).append(data)
    result={}
    for mode,runs in groups.items():
        rows=[r['systems'] for r in runs]
        factor=[sum(x['factor_seconds'] for x in row) for row in rows]
        solve=[sum(x['two_rhs_seconds'] for x in row) for row in rows]
        actual=lambda r:[s.get('actual_q_fnv1a64',r.get('actual_q_fnv1a64')) for s in r['systems']]
        stats=['lnz','unz','fill','flops','off_diagonal_pivots']
        result[mode]=dict(samples=len(runs),factor_seconds=factor,two_rhs_seconds=solve,
            median_factor_seconds=statistics.median(factor),
            median_service_seconds=statistics.median([a+b+r['analysis_seconds'] for a,b,r in zip(factor,solve,runs)]),
            factor_range_over_mean_percent=100*(max(factor)-min(factor))/statistics.mean(factor),
            mean_statistics={key:(statistics.mean(x[key] for x in rows[0] if x[key]>=0)
                if any(x[key]>=0 for x in rows[0]) else None) for key in stats if key in rows[0][0]},
            actual_q_hashes=actual(runs[0]),
            repeated_factors_and_permutations_identical=all(
                [{k:x.get(k) for k in stats} for x in row]==[{k:x.get(k) for k in stats} for x in rows[0]] and actual(r)==actual(runs[0])
                for row,r in zip(rows,runs)),
            all_match_given_q=all(x['matches_input_q'] for x in rows[0]) if 'matches_input_q' in rows[0][0] and rows[0][0]['matches_input_q'] is not None else None,
            worst_raw_componentwise_backward_error=max(x['quality']['raw']['componentwise_backward_error'] for row in rows for x in row),
            worst_raw_normwise_backward_error=max(x['quality']['raw']['normwise_backward_error'] for row in rows for x in row))
    comparisons=[]
    pairs=[('eigen_colamd','eigen_amd'),('default','unsym'),('default','noscale'),('default','pivot1'),
           ('given_colamd','given_amd'),('given_colamd','given_colamd_noscale'),
           ('given_colamd','given_colamd_unsym'),('given_colamd','given_colamd_pivot1'),
           ('eigen_colamd','given_colamd_noscale'),('eigen_colamd','given_colamd_pivot1')]
    for a,b in pairs:
        if a not in result or b not in result:continue
        x,y=result[a],result[b]
        comparisons.append(dict(baseline=a,candidate=b,same_actual_q=x['actual_q_hashes']==y['actual_q_hashes'],
            factor_saving_percent=100*(1-y['median_factor_seconds']/x['median_factor_seconds']),
            fill_change_percent=100*(y['mean_statistics']['fill']/x['mean_statistics']['fill']-1)))
    return dict(status=summary['status'],variants=result,comparisons=comparisons,
        scope='Frozen historical systems. Median of repeated numeric factor totals; quality and diagnostics outside timing. Not curve qualification or exclusive BLAS attribution.')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path);a=p.parse_args()
    result=summarize(a.root);write(a.root/'analysis.json',result);print(json.dumps(result,indent=2))
