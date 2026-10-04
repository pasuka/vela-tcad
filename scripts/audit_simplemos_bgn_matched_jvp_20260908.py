"""Repeat the frozen block-Jv diagnostic at R1b matched-ni states."""
import argparse
from pathlib import Path
import numpy as np
import simplemos_bgn_intrinsic_control_20260908 as b

w=b.w
a,d,V,v=b.a,b.d,b.V,b.v
LOCAL=b.LOCAL/'jvp'
OUT=b.OUT/'jvp'
STEPS=(1e-4,3e-5,1e-5,3e-6)


def prepare():
    a.verify(b.OUT/'comparison_evidence.json')
    jobs=[]
    for c in a.read(b.OUT/'contract.json')['cases']:
        root=b.LOCAL/'vela'/c['case']/'vela/vg_040'
        src=next(p for p in sorted(root.glob('attempt_*')) if a.read(p/'result.json')['qualified'])
        geo,mask=V.m.previous.prior.support(c)
        dest=LOCAL/c['case']
        cfg=a.read(src/'config.json')
        cfg.update(simulation_type='newton_jvp_probe',state_file=str(src/'state.csv'),output_csv=str(dest/'jvp.csv'),
            directions=[dict(name=f'{mode}_{h:.0e}',mode=mode,amplitude_V=h,node_ids=np.where(mask)[0][::3].tolist(),exclude_contacts=True)
                        for mode in ('psi','phin','phip') for h in STEPS])
        cfg.pop('output_state_file',None)
        a.write(dest/'config.json',cfg)
        jobs.append(dict(case=c['case'],device=c['device'],vd=c['vd'],config=str(dest/'config.json'),state=str(src/'state.csv')))
    a.write(OUT/'contract.json',dict(jobs=jobs,steps_V=STEPS,vg=.8,
        purpose='Read-only derivative diagnostic at matched OldSlotboom intrinsic concentration. Does not change DC acceptance or repair the Jacobian.',
        gate='Same true per-block relative threshold 1e-4 as R1; an extra 3e-6 V amplitude tests the error plateau. Cross-carrier source blocks remain separately labeled.',
        hypotheses='BGN uses the variable-ni SG path; production stable psi-derivative replacement is currently restricted to no-BGN/equal-ni. A failed Jv is not alone proof of the failure mechanism.'))
    d.matrix.freeze(OUT/'freeze.json',[Path(__file__).resolve(),b.OUT/'comparison_evidence.json',OUT/'contract.json',w.RUNNER]+
                    [Path(job[k]) for job in jobs for k in ('config','state')])


def run():
    a.verify(OUT/'freeze.json')
    rows=[]
    for job in a.read(OUT/'contract.json')['jobs']:
        path=Path(job['config'])
        status=V.execute(path,w.RUNNER,V.environment())
        assert status['exit_code']==0
        for r in a.rows(path.parent/'jvp.csv'):
            for block in ('psi','phin','phip'):
                ana,fd=float(r[f'analytic_{block}_norm']),float(r[f'finite_difference_{block}_norm'])
                absolute=float(r[f'{block}_relative_error'])*max(1.,fd)
                relative=absolute/max(ana,fd,1e-300)
                weak=(r['mode'],block) in (('phin','phip'),('phip','phin'))
                rows.append(dict(case=job['case'],device=job['device'],vd=job['vd'],vg=.8,input_block=r['mode'],
                    output_block=block,amplitude_V=float(r['amplitude_V']),analytic_norm=ana,fd_norm=fd,
                    absolute_error=absolute,true_relative=relative,weak_carrier_cross_block=weak,qualified=relative<=1e-4))
    v.csv_union(OUT/'blocks.csv',rows)
    failed=[r for r in rows if not r['weak_carrier_cross_block'] and r['amplitude_V']<1e-4 and not r['qualified']]
    a.write(OUT/'summary.json',dict(states=4,directions=48,block_checks=len(rows),
        smaller_amplitude_strong_checks=84,smaller_amplitude_strong_failures=len(failed),
        derivative_gate_passed=not failed,max_true_relative=max(r['true_relative'] for r in rows)))
    d.matrix.freeze(OUT/'evidence.json',[OUT/'freeze.json',OUT/'blocks.csv',OUT/'summary.json']+[p for p in LOCAL.rglob('*') if p.is_file()])
    print(a.read(OUT/'summary.json'),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('prepare','run'))
    globals()[parser.parse_args().action]()
