"""Run the first double initialization only after its own response/column gates."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import validate_simplemos_enormal_candidate_20260912 as q
import simplemos_enormal_material_scale_20260912 as n

def main():
    a,d=q.a,q.d
    a.verify(q.O/'candidate_freeze.json');a.verify(q.O/'fixed_evidence.json');a.verify(n.O/'pilot_evidence.json')
    assert a.read(n.O/'pilot_summary.json')['all_qualified']
    cc=a.read(q.O/'inputs.json')['cases'][0]
    assert cc['case']==a.read(n.O/'native_contract.json')['pilot_case']
    dest=q.L/'real_columns'/cc['key'];a.verify(q.O/'smoke_evidence.json')
    status=a.read(dest/'jvp.status.json');assert status['exit_code']==0
    rows=a.rows(dest/'rows.csv');assert len(rows)==26640
    import math
    errors=[]
    for r in rows:
        an=float(r['analytic_derivative']);fd=float(r['finite_difference_derivative'])
        assert math.isfinite(an) and math.isfinite(fd)
        scale=max(abs(an),abs(fd));errors.append(abs(an-fd)/scale if scale else 0.)
    assert max(errors)<=1e-4
    freeze=q.O/'pilot_dc_prerequisites.json'
    if not freeze.exists():d.matrix.freeze(freeze,[Path(__file__).resolve(),q.O/'candidate_freeze.json',q.O/'fixed_evidence.json',n.O/'pilot_evidence.json',n.O/'pilot_summary.json']+list(dest.iterdir()))
    a.verify(freeze);q.configure()
    def run(arm):return q.v.solve_target(cc,cc['index'],arm,Path(cc['phumob_seed' if arm=='phumob_seed' else 'native_seed']))
    with ThreadPoolExecutor(max_workers=2) as pool:attempts=[r for group in pool.map(run,('phumob_seed','native_seed')) for r in group]
    q.v.csv_union(q.O/'pilot_dc_attempts.csv',attempts)
    d.matrix.freeze(q.O/'pilot_dc_evidence.json',[freeze,q.O/'pilot_dc_attempts.csv']+[f for f in (q.L/'vela'/cc['case']).rglob('*') if f.is_file()])
    print(attempts,flush=True)

if __name__=='__main__':main()
