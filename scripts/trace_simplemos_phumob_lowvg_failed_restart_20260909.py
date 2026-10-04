"""Three-iteration production trace; retains refinement, gates and failures."""
from pathlib import Path
import calibrate_simplemos_phumob_local_source_20260909 as s

LOCAL=s.LOCAL/'failed_trace';OUT=s.OUT/'failed_trace'

def run():
    case='m65_n23_vd_1p000000_endpoint';root=s.prior.LOCAL/'candidate/dc/phumob'/case
    jobs=[];files=[Path(__file__).resolve(),s.prior.q.run.RUNNER,s.REPO/'src/solver/NewtonSolver.cpp']
    for label,relative in (('native_first','native/vg_010/attempt_0'),('native_reload','native/vg_010/attempt_1')):
        base=root/relative;dest=LOCAL/label;cfg=s.a.read(base/'config.json')
        cfg.update(state_file=str(base/'state.csv'),output_state_file=str(dest/'state.csv'))
        cfg['solver']['max_iter']=3;cfg['solver']['carrier_row_convergence']['min_newton_max_iter']=3
        cfg['solver']['local_update_diagnostics']=dict(enabled=True,nodes=[967,983,1089],csv_file=str(dest/'updates.csv'),first_iterations=3,every_iterations=1)
        assert cfg['solver']['linear_refinement_iterations']==4 and cfg['solver']['carrier_row_convergence']['eps_row']==1e-6
        s.a.write(dest/'config.json',cfg);s.V.post_config(cfg,dest)
        jobs.append(dict(label=label,base=str(base),dest=str(dest)));files += [base/'state.csv',base/'config.json']+list(dest.glob('*.json'))
    s.a.write(OUT/'contract.json',dict(jobs=jobs,scope='Budget-limited trajectory diagnostic from previous failure outputs; not a third qualification retry. Same production binary, scaling, caps, four linear refinements, line search and acceptance eps. max_iter and carrier minimum iteration budget capped at 3.',prior_step_probe_limitation='evaluateStep uses a direct default LinearSolver, omits solve() row equilibration and linear refinements, and evaluates one full trial without line search. Its direction is not the exact iterative solver direction. This trace supplies actual production updates.'))
    s.d.matrix.freeze(OUT/'freeze.json',files+[OUT/'contract.json'])
    results=[]
    for job in jobs:
        dest=Path(job['dest']);status=s.V.execute(dest/'config.json',s.prior.q.run.RUNNER,s.V.environment())
        for name in ('all_row','acceptance_edges'):s.V.execute(dest/(name+'.json'),s.prior.q.run.RUNNER,s.V.environment())
        result=s.prior.q.run.w.old.prior.old.qualify(dict(device='n23'),dest)
        results.append(dict(label=job['label'],**result));print(results[-1],flush=True)
    s.a.write_csv(OUT/'summary.csv',results);s.d.matrix.freeze(OUT/'evidence.json',[OUT/'freeze.json',OUT/'summary.csv']+[p for p in LOCAL.rglob('*') if p.is_file()])

if __name__=='__main__':run()
