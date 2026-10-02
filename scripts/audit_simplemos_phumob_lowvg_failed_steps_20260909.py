"""Read-only raw/clipped Newton-step checks on both failed states and control."""
from decimal import Decimal, localcontext
from pathlib import Path
import calibrate_simplemos_phumob_local_source_20260909 as s

LOCAL=s.LOCAL/'failed_steps';OUT=s.OUT/'failed_steps'

def run():
    case='m65_n23_vd_1p000000_endpoint';root=s.prior.LOCAL/'candidate/dc/phumob'/case
    jobs=[];files=[Path(__file__).resolve(),s.prior.q.run.RUNNER,s.REPO/'src/solver/NewtonSolver.cpp']
    for label,relative in (('native_first','native/vg_010/attempt_0'),('native_reload','native/vg_010/attempt_1'),('qualified_control','vela/vg_010/attempt_0')):
        base=root/relative;dest=LOCAL/label;cfg=s.a.read(base/'config.json')
        cfg.update(simulation_type='newton_step_probe',state_file=str(base/'state.csv'),output_csv=str(dest/'step.csv'));cfg.pop('output_state_file',None)
        s.a.write(dest/'config.json',cfg);jobs.append(dict(label=label,base=str(base),dest=str(dest)));files += [base/'state.csv',base/'all_row.csv',base/'config.status.json',dest/'config.json']
    s.a.write(OUT/'contract.json',dict(jobs=jobs,scope='Evaluate one current-production Newton direction and clipped full trial on frozen final states. Does not replay line search or accept a new state. Original all-row gate 1e-6 retained.',norm='Decimal100 sum(trial residual squared - baseline residual squared) using exactly the exported double residual values. Evaluates norm comparison, not residual-kernel precision.'))
    s.d.matrix.freeze(OUT/'freeze.json',files+[OUT/'contract.json'])
    results=[];hotspots=[]
    for job in jobs:
        dest=Path(job['dest']);base=Path(job['base']);status=s.V.execute(dest/'config.json',s.prior.q.run.RUNNER,s.V.environment());assert status['exit_code']==0,status
        rows=s.a.rows(dest/'step.csv');terms=s.d.ordered(base/'all_row.csv',len(rows));old=s.a.read(base/'config.status.json')
        with localcontext() as ctx:
            ctx.prec=100
            differences={block:sum((Decimal(r['trial_'+block+'_residual'])**2-Decimal(r[block+'_residual'])**2 for r in rows),Decimal(0)) for block in ('psi','phin','phip')}
            total=sum(differences.values(),Decimal(0))
        for v in old['carrier_row_convergence'].get('violations',[]):
            i=v['node_id'];r=rows[i];car=v['carrier'];block='phin' if car=='electron' else 'phip';scale=max(float(terms[i][car+'_flux_abs_sum']),abs(float(terms[i][car+'_recombination'])),abs(float(terms[i][car+'_impact'])))
            hotspots.append(dict(label=job['label'],node=i,carrier=car,original_ratio=v['ratio'],reloaded_residual=float(r[block+'_residual']),trial_residual=float(r['trial_'+block+'_residual']),trial_over_old_scale=abs(float(r['trial_'+block+'_residual']))/scale,delta_psi_V=r['delta_psi_V'],delta_qf_V=r['delta_'+block+'_V'],psi_update_absorbed=float(r['trial_psi'])==float(r['psi']),trial_is_accepted=False))
        results.append(dict(label=job['label'],original_qualified=old['converged'],original_max_row=old['carrier_row_convergence']['max_ratio'],original_worst_node=old['carrier_row_convergence']['max_ratio_node'],original_worst_carrier=old['carrier_row_convergence']['max_ratio_carrier'],raw_step_norm=status['raw_step_norm'],clipped_step_norm=status['step_norm'],clipping_changed=status['raw_step_norm']!=status['step_norm'],double_norm_decreases=status['trial_block_residuals']['combined']<status['block_residuals']['combined'],exact_exported_norm_decreases=total<0,exact_exported_norm_squared_difference=str(total),**{k+'_norm_squared_difference':str(v) for k,v in differences.items()}))
        print(results[-1],flush=True)
    s.a.write_csv(OUT/'summary.csv',results);s.a.write_csv(OUT/'violations.csv',hotspots)
    s.d.matrix.freeze(OUT/'evidence.json',[OUT/'freeze.json',OUT/'summary.csv',OUT/'violations.csv']+[p for p in LOCAL.rglob('*') if p.is_file()])

if __name__=='__main__':run()
