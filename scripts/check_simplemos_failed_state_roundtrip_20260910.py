"""Zero-iteration state roundtrip, diagnostic only and no acceptance credit."""
import copy
from pathlib import Path
import validate_simplemos_production_weighted_20260910 as v

def main():
    a,d=v.a,v.d;local=v.REPO/'build-release/pp_zero';out=v.OUT/'failed_roundtrip'
    a.verify(v.OUT/'failed_restart_isolation/evidence.json')
    base=v.LOCAL/'dc/phumob/m65_n23_vd_1p000000_endpoint/native/vg_000/attempt_1'
    cfg=a.read(base/'config.json');job=dict(case='m65_n23_vd_1p000000_endpoint',device='n23',vd=1.,vg=0.)
    geo,mask=v.run.V.m.previous.prior.support(job);seed=base/'state.csv';rows=[];files=[Path(__file__).resolve(),v.OUT/'failed_restart_isolation/evidence.json',seed]
    for k in range(2):
        dest=local/str(k);deck=copy.deepcopy(cfg);deck.update(state_file=str(seed),output_state_file=str(dest/'state.csv'))
        deck['solver']['max_iter']=0;deck['solver']['carrier_row_convergence']['mode']='report'
        a.write(dest/'config.json',deck);s=v.run.V.execute(dest/'config.json',v.RUNNER,v.run.V.environment());assert s['iterations']==0,s
        delta=v.run.V.m.previous.prior.delta_states(d.ordered(seed,geo.count),d.ordered(dest/'state.csv',geo.count),mask)
        row=dict(reload=k+1,initial_residual=s['initial_residual'],final_residual=s['final_residual'],**delta)
        row['psi_changed_nodes']=sum(float(x['psi'])!=float(y['psi']) for x,y in zip(d.ordered(seed,geo.count),d.ordered(dest/'state.csv',geo.count)))
        rows.append(row);seed=dest/'state.csv';files += [p for p in dest.iterdir() if p.is_file()]
    old=a.read(base/'config.status.json')
    summary=dict(prior_final_residual=old['final_residual'],first_reload_initial_residual=rows[0]['initial_residual'],
        norm_ratio=rows[0]['initial_residual']/old['final_residual'],zero_Newton_iterations=True,cohort_credit=False,
        interpretation='A saved-state roundtrip changes the residual before any Newton update. This diagnoses numerical representation/repacking; it does not establish a complete remedy or justify additional cohort reloads.')
    a.write_csv(out/'roundtrip.csv',rows);a.write(out/'summary.json',summary)
    d.matrix.freeze(out/'evidence.json',files+[out/'roundtrip.csv',out/'summary.json']);print(rows,summary,flush=True)

if __name__=='__main__':main()
