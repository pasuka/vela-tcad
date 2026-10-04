"""Isolate the plain PhuMob Jacobian repair with frozen inputs and old gates."""
import argparse,copy,math
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
import audit_simplemos_phumob_assembled_cross_20260909 as cross
import analyze_simplemos_phumob_assembled_cross_20260909 as analysis

p=cross.p;q=p.qualified;a,d=p.a,p.d
LOCAL=p.REPO/'build-release/phumob_fix_20260909'
OUT=p.REPO/'reference_tcad/simplemos_sentaurus2022/phumob_fix_20260909'
BASELINE=LOCAL/'baseline/vela_example_runner.exe'
RUNNER=q.run.RUNNER
OLD_CROSS=cross.OUT


def prepare():
    a.verify(OUT/'baseline.json')
    jobs=[];files=[Path(__file__).resolve(),OUT/'baseline.json']
    for old in a.read(q.OUT/'validation_contract.json')['jobs']:
        if old['model']!='old_slotboom' or old['index']!=40:continue
        job=dict(old);job['model']='phumob'
        cfg=a.read(Path(old['original_config']))
        cfg['solver']['mobility'].update(model='phumob',edge_averaging='legacy')
        dest=LOCAL/'inputs'/old['case']/old['arm'];path=dest/'config.json'
        a.write(path,cfg);job['original_config']=str(path);jobs.append(job)
        files += [path,Path(job['seed'])]+[Path(cfg[k]) for k in ('mesh_file','materials_file','node_doping_file')]
    assert len(jobs)==8
    a.write(OUT/'contract.json',dict(jobs=jobs,gates=q.run.w.GATES,scope='Vg=.8, 2 NWell x 2 Vd x two inherited independent initializations. Before/after same PhuMob legacy constitutive and same SRH/BGN; only Jacobian changes.',
        recovery='At most one same-bias saved-state reload; all attempts retained. A failed point is not used for qualified final-state comparison.',
        caveat='This is a derivative-isolation self-consistent comparison, not native element-box PhuMob restoration. Native current is context only while the averaging definition differs.'))
    d.matrix.freeze(OUT/'input_freeze.json',files+[OUT/'contract.json'])


def solve(arm):
    a.verify(OUT/'input_freeze.json')
    runner=BASELINE if arm=='before' else RUNNER
    if arm=='after':
        a.verify(OUT/'cross_analysis/evidence.json')
        summary=a.read(OUT/'cross_analysis/summary.json');assert summary['derivative_failures']==0 and summary['reference_units_and_state_qualified']
        assert a.read(OUT/'blocks_summary.json')['gated_failures']==0
    q.run.LOCAL=LOCAL/arm;q.run.OUT=OUT/arm;q.run.RUNNER=runner
    jobs=a.read(OUT/'contract.json')['jobs']
    a.write(q.run.OUT/'validation_contract.json',dict(jobs=jobs))
    d.matrix.freeze(q.run.OUT/'validation_freeze.json',[OUT/'input_freeze.json',runner,OUT/'baseline.json'])
    with ThreadPoolExecutor(max_workers=2) as pool:
        rows=[r for group in pool.map(q.run.solve_job,jobs) for r in group]
    q.run.v.csv_union(q.run.OUT/'attempts.csv',rows)
    d.matrix.freeze(q.run.OUT/'dc_evidence.json',[q.run.OUT/'validation_freeze.json',q.run.OUT/'attempts.csv']+[x for x in (q.run.LOCAL/'dc').rglob('*') if x.is_file()])


def cross_check():
    a.verify(OUT/'input_freeze.json')
    contract=a.read(OLD_CROSS/'contract.json')
    jobs=[];files=[OUT/'input_freeze.json',RUNNER]
    for old in contract['jobs']:
        job=dict(old);dest=LOCAL/'cross'/old['case'];job['dir']=str(dest);jobs.append(job)
        for probe in ('edges','jvp'):
            original=Path(old['dir'])/(probe+'.json');cfg=a.read(original)
            cfg['output_csv']=str(dest/(probe+'.csv'))
            if probe=='jvp':cfg['row_output_csv']=str(dest/'rows.csv')
            a.write(dest/(probe+'.json'),cfg);files += [original,dest/(probe+'.json'),Path(cfg['state_file'])]
    cross.LOCAL=LOCAL/'cross';cross.OUT=OUT/'cross';contract['jobs']=jobs
    a.write(cross.OUT/'contract.json',contract)
    files += [cross.OUT/'contract.json']+[p.REPO/n for n in ('src/equation/CoupledDDAssembler.cpp','src/physics/MobilityModel.cpp','include/vela/physics/MobilityModel.h','include/vela/equation/AssemblerUtils.h')]
    d.matrix.freeze(cross.OUT/'freeze.json',files)
    cross.execute()
    analysis.LOCAL=cross.LOCAL;analysis.OUT=OUT/'cross_analysis'
    d.matrix.freeze(analysis.OUT/'analysis_freeze.json',[Path(analysis.__file__),cross.OUT/'raw_evidence.json'])
    analysis.analyze()


def blocks():
    a.verify(OUT/'cross_analysis/evidence.json')
    jobs=[];files=[]
    for old in a.read(OLD_CROSS/'contract.json')['jobs']:
        src=Path(old['dir'])/'edges.json';cfg=a.read(src)
        _,mask=q.run.V.m.previous.prior.support(old)
        dest=LOCAL/'blocks'/old['case'];cfg.update(simulation_type='newton_jvp_probe',output_csv=str(dest/'jvp.csv'))
        cfg['directions']=[dict(name=f'{mode}_{step:.0e}',mode=mode,amplitude_V=step,node_ids=np.where(mask)[0][::3].tolist(),exclude_contacts=True)
            for mode in ('psi','phin','phip') for step in (1e-4,3e-5,1e-5,3e-6)]
        a.write(dest/'config.json',cfg);jobs.append(dict(case=old['case'],path=str(dest/'config.json')));files += [src,dest/'config.json',Path(cfg['state_file'])]
    a.write(OUT/'blocks_contract.json',dict(jobs=jobs,gate=1e-4,gated_steps='Smaller three amplitudes, excluding weak cross blocks already compared to independent HP references. All outputs retained.'))
    d.matrix.freeze(OUT/'blocks_freeze.json',files+[OUT/'blocks_contract.json',RUNNER,Path(__file__).resolve()])
    def one(job):
        path=Path(job['path']);status=q.run.V.execute(path,RUNNER,q.run.V.environment());assert status['exit_code']==0,status
        result=[]
        for row in a.rows(path.parent/'jvp.csv'):
            for block in ('psi','phin','phip'):
                ana=float(row[f'analytic_{block}_norm']);fd=float(row[f'finite_difference_{block}_norm'])
                error=float(row[f'{block}_relative_error'])*max(1.,fd);relative=error/max(ana,fd,1e-300)
                weak=(row['mode'],block) in (('phin','phip'),('phip','phin'))
                result.append(dict(case=job['case'],input_block=row['mode'],output_block=block,step_V=row['amplitude_V'],true_relative=relative,qualified=relative<=1e-4,gated=not weak and float(row['amplitude_V'])<1e-4))
        return result
    with ThreadPoolExecutor(max_workers=2) as pool:rows=[r for group in pool.map(one,jobs) for r in group]
    a.write_csv(OUT/'blocks.csv',rows);gated=[r for r in rows if r['gated']]
    summary=dict(checks=len(rows),gated_checks=len(gated),gated_failures=sum(not r['qualified'] for r in gated),gated_max_relative=max(r['true_relative'] for r in gated))
    a.write(OUT/'blocks_summary.json',summary)
    d.matrix.freeze(OUT/'blocks_evidence.json',[OUT/'blocks_freeze.json',OUT/'blocks.csv',OUT/'blocks_summary.json']+[x for x in (LOCAL/'blocks').rglob('*') if x.is_file()])
    print(summary,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('prepare','before','cross_check','blocks','after'))
    action=parser.parse_args().action
    solve(action) if action in ('before','after') else globals()[action]()
