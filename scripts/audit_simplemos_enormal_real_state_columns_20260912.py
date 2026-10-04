"""True per-entry real-mesh column checks supplement the complete small-mesh tests."""
import math,time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import validate_simplemos_enormal_candidate_20260912 as q

a,d,L,O=q.a,q.d,q.L,q.O

def main():
    a.verify(O/'candidate_freeze.json');a.verify(O/'smoke_evidence.json')
    cases=a.read(O/'inputs.json')['cases']
    def run(cc):
        cfg=q.config_case(cc);count=len(a.read(Path(cfg['mesh_file']))['nodes'])
        dest=L/'real_columns'/cc['key']
        cfg.update(simulation_type='newton_jvp_probe',state_file=cc['native_seed'],output_csv=str(dest/'jvp.csv'),row_output_csv=str(dest/'rows.csv'),
            directions=[dict(name=f'{mode}_{h:.0e}',mode=mode,amplitude_V=h,node_ids=[1000],exclude_contacts=False) for mode in ('psi','phin','phip') for h in (1e-6,1e-20)],
            sample_rows=[dict(block=b,node_id=i) for b in ('psi','phin','phip') for i in range(count)])
        a.write(dest/'jvp.json',cfg);s=q.V.execute(dest/'jvp.json',q.RUNNER,q.V.environment());assert s['exit_code']==0,s
        groups={};failures=[]
        for row in a.rows(dest/'rows.csv'):
            an=float(row['analytic_derivative']);fd=float(row['finite_difference_derivative'])
            finite=math.isfinite(an) and math.isfinite(fd);scale=max(abs(an),abs(fd));rel=abs(an-fd)/scale if scale else 0.
            key=(row['direction'],row['row_block']);g=groups.setdefault(key,dict(case=cc['key'],direction=key[0],row_block=key[1],rows=0,nonzero_entries=0,max_relative=0.,qualified=True))
            g['rows']+=1;g['nonzero_entries']+=scale>0;g['max_relative']=max(g['max_relative'],rel)
            good=finite and rel<=1e-4;g['qualified'] &= good
            if not good:failures.append(dict(case=cc['key'],**row,true_relative_error=rel))
        result=dict(case=cc['key'],elapsed_seconds=s['elapsed_seconds'],checks=sum(r['rows'] for r in groups.values()),failed=len(failures),qualified=not failures)
        print(result,flush=True)
        return list(groups.values()),failures,result
    with ThreadPoolExecutor(max_workers=2) as pool:results=list(pool.map(run,cases))
    rows=[r for g,_,_ in results for r in g];failures=[r for _,f,_ in results for r in f]
    a.write_csv(O/'real_column_groups.csv',rows)
    a.write(O/'real_column_failures.json',failures)
    summary=dict(cases=[r for _,_,r in results],qualified=all(r['qualified'] for _,_,r in results),scope='All residual rows for the three columns of node 1000, at two amplitudes, across eight real states. This does not certify every column of each real-state matrix; all columns and weak nonlocal entries are additionally checked on the controlled synthetic mesh.',gate=1e-4,excluded_weak_entries=0)
    a.write(O/'real_column_summary.json',summary)
    d.matrix.freeze(O/'real_column_evidence.json',[Path(__file__).resolve(),O/'candidate_freeze.json',O/'real_column_groups.csv',O/'real_column_failures.json',O/'real_column_summary.json']+[f for f in (L/'real_columns').rglob('*') if f.is_file()])
    assert summary['qualified'],summary

if __name__=='__main__':main()
