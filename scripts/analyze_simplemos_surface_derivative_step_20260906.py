"""Independent classification of frozen local transport derivative step controls."""
from pathlib import Path
import audit_simplemos_surface_derivative_step_20260906 as s

a=s.a;m=s.m


def main():
    a.verify(s.OUT/'freeze.json')
    rows=[];invariance=[];files=[Path(__file__).resolve(),s.OUT/'freeze.json']
    for job in a.read(s.OUT/'contract.json')['jobs']:
        path=Path(job['config']);status=path.with_suffix('.status.json')
        assert a.read(status)['exit_code']==0
        raw=a.rows(path.parent/'columns.csv');files += [status,path.parent/'columns.csv',path.parent/'columns.csv.residual.csv']
        baseline=m.LOCAL/job['case']/'baseline_on/columns.csv.residual.csv';files.append(baseline)
        same=a.sha(baseline)==a.sha(path.parent/'columns.csv.residual.csv');assert same
        invariance.append(dict(case=job['case'],transport_step=job['step'],residual_identical=same))
        primary=[r for r in raw if abs(float(r['step_V'])/1e-6-1)<1e-10]
        fine={(r['column_node'],r['column_block'],r['row_block']):r for r in raw if abs(float(r['step_V'])/1e-7-1)<1e-10}
        assert len(primary)==162
        for cb in range(3):
            for rb in range(3):
                group=[r for r in primary if int(r['column_block'])==cb and int(r['row_block'])==rb]
                peak=max(float(r['fd_norm']) for r in group)
                for r in group:
                    f=fine[r['column_node'],r['column_block'],r['row_block']];signal=float(r['fd_norm'])
                    active=signal>max(1e-10*peak,1e-300)
                    stable=active and float(f['fd_step_difference_norm'])<=1e-4*signal
                    error=float(r['diff_norm'])/max(signal,1e-300)
                    rows.append(dict(case=job['case'],transport_step=job['step'],column_node=r['column_node'],column_block=cb,row_block=rb,
                        active=active,fd_stable=stable,relative_error=error,fine_relative_error=float(f['diff_norm'])/max(float(f['fd_norm']),1e-300),
                        passed=stable and error<=1e-3,max_diff_node=r['max_diff_node'],outside_pattern_relative=float(r['outside_pattern_norm'])/max(signal,1e-300)))
    a.write_csv(s.OUT/'columns.csv',rows);a.write_csv(s.OUT/'invariance.csv',invariance)
    fails=[r for r in a.rows(m.OUT/'columns.csv') if r['variant']=='baseline_on' and r['fd_stable']=='True' and r['passed']=='False']
    controls=[]
    for f in fails:
        matches=[r for r in rows if (r['case'],str(r['column_node']),str(r['column_block']),str(r['row_block']))==(f['case'],f['column_node'],f['column_block'],f['row_block'])]
        assert len(matches)==2
        for r in matches:controls.append(dict(**r,baseline_relative_error=float(f['relative_error'])))
    a.write_csv(s.OUT/'defect_controls.csv',controls)
    a.write(s.OUT/'result.json',dict(input_hashes={a.rel(p):a.sha(p) for p in sorted(set(files))},
        stable_failures_remaining=sum(r['fd_stable'] and not r['passed'] for r in rows),
        previously_failed_controls_passed=sum(r['passed'] for r in controls),previously_failed_controls=len(controls),
        residual_invariant=True,production_promoted=False,
        interpretation='Local transport FD step control, not a blanket certificate for all Jacobian entries or a production fix.'))
    print(controls,flush=True)


if __name__=='__main__':main()
