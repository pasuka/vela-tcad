"""Validate the BGN psi-derivative repair with frozen state and DC controls."""
import argparse
import copy
import math
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
import simplemos_bgn_intrinsic_control_20260908 as b
import audit_simplemos_bgn_matched_jvp_20260908 as old_jvp
import audit_simplemos_bgn_psi_precision_20260908 as precision

w,v,V,a,d=b.w,b.v,b.V,b.a,b.d
LOCAL,OUT=precision.LOCAL,precision.OUT
RUNNER=w.RUNNER
STEPS=(1e-4,3e-5,1e-5,3e-6)


def prepare():
    a.verify(OUT/'preimage_freeze.json');a.verify(OUT/'precision_evidence.json')
    a.verify(b.OUT/'comparison_evidence.json');a.verify(old_jvp.OUT/'evidence.json')
    jobs=[];files=[Path(__file__).resolve(),RUNNER,w.REPO/'build-release/libvela_core.a',OUT/'precision_evidence.json',b.OUT/'comparison_evidence.json']
    for c in a.read(b.OUT/'contract.json')['cases']:
        for model,indices,root in (('old_slotboom',w.n.INDICES,b.LOCAL/'vela'/c['case']),
                                   ('no_bgn',(40,50),w.LOCAL/'vela/no_bgn'/c['case'])):
            for index in indices:
                for arm in ('vela','native'):
                    src=root/arm/f'vg_{index:03d}'/'attempt_0/config.json'
                    cfg=a.read(src)
                    seed=Path(cfg['state_file'])
                    jobs.append(dict(case=c['case'],device=c['device'],vd=c['vd'],vg=w.n.prior.GRID[index],index=index,
                                     model=model,arm=arm,original_config=str(src),seed=str(seed)))
                    files += [src,seed,Path(cfg['materials_file']),Path(cfg['mesh_file']),Path(cfg['node_doping_file'])]
    sources=[w.REPO/p for p in ('src/equation/CoupledDDAssembler.cpp','include/vela/discretization/StableSGDerivative.h',
                               'tests/test_production_numerics.cpp','tests/test_sg_flux.cpp')]
    a.write(OUT/'validation_contract.json',dict(jobs=jobs,points=24,initial_attempts=48,gates=w.GATES,
        physics_changed=False,acceptance_changed=False,
        change='Fixed-mobility BGN psi partials use the already stable edge flux and its effective Bernoulli logarithmic derivative; no residual or current formula changes.',
        recovery='At most one same-bias saved-state reload after a failed attempt; retain both attempts.',
        jvp=dict(steps_V=STEPS,relative_gate=1e-4,
                 same_state='Four original matched-ni Vg=.8 states and exact original directions; compare with frozen pre-repair Jv ledger.',
                 solved_state='All 24 post-repair targets, every third free Si node, psi/phin/phip separately; use both smaller steps for the gate.',
                 weak_cross='Label separately; full-residual differences do not replace isolated weak source-block completeness tests.'),
        further_models='PhuMob, Enormal and high-field saturation remain disabled until this repair and their independent cell/edge constitutive calibration pass.'))
    d.matrix.freeze(OUT/'validation_freeze.json',files+sources+[OUT/'validation_contract.json',Path(w.__file__),Path(b.__file__),Path(precision.__file__)])
    print('Frozen 16 BGN + 8 no-BGN points x two initializations.',flush=True)


def attempt(job,seed,number):
    dest=LOCAL/'dc'/job['model']/job['case']/job['arm']/f"vg_{job['index']:03d}"/f'attempt_{number}'
    cfg=a.read(Path(job['original_config']))
    cfg.update(state_file=str(seed),output_state_file=str(dest/'state.csv'))
    if not (dest/'input_freeze.json').exists():
        a.write(dest/'config.json',cfg);V.post_config(cfg,dest)
        fun=copy.deepcopy(a.read(dest/'acceptance_edges.json'))
        fun.update(simulation_type='terminal_current_functional_probe',contact='drain',
                   residual_output_csv=str(dest/'port_residual.csv'),contact_edge_output_csv=str(dest/'port_edges.csv'))
        fun.pop('output_csv',None);a.write(dest/'functional.json',fun)
        d.matrix.freeze(dest/'input_freeze.json',[Path(seed),dest/'config.json',dest/'all_row.json',dest/'acceptance_edges.json',dest/'functional.json',OUT/'validation_freeze.json'])
    else:assert a.read(dest/'config.json')==cfg
    a.verify(dest/'input_freeze.json')
    status=V.execute(dest/'config.json',RUNNER,V.environment())
    row=dict(**job,attempt=number,dest=str(dest),qualified=False,exit_code=status['exit_code'],
             failure=status.get('failure_reason',''),elapsed_seconds=status['elapsed_seconds'])
    if (dest/'state.csv').exists():
        for name in ('all_row','acceptance_edges','functional'):V.execute(dest/(name+'.json'),RUNNER,V.environment())
        row.update(w.old.prior.old.qualify(job,dest))
        port=a.read(dest/'functional.status.json')
        error=abs(port['current_A_per_um']/port['contact_current_extractor_A_per_um']-1)
        row.update(port_relative=error,qualified=row['qualified'] and port['exit_code']==0 and error<=w.GATES['port_relative'])
    v.write_same(dest/'result.json',row)
    print(job['model'],job['device'],job['vd'],job['vg'],job['arm'],number,row['qualified'],'row',row.get('max_row_ratio'),flush=True)
    return row


def solve_job(job):
    rows=[attempt(job,Path(job['seed']),0)]
    if not rows[-1]['qualified'] and (Path(rows[-1]['dest'])/'state.csv').exists():
        rows.append(attempt(job,Path(rows[-1]['dest'])/'state.csv',1))
    return rows


def solve():
    a.verify(OUT/'validation_freeze.json')
    with ThreadPoolExecutor(max_workers=2) as pool:
        rows=[r for block in pool.map(solve_job,a.read(OUT/'validation_contract.json')['jobs']) for r in block]
    v.csv_union(OUT/'attempts.csv',rows)
    d.matrix.freeze(OUT/'dc_evidence.json',[OUT/'validation_freeze.json',OUT/'attempts.csv']+[p for p in (LOCAL/'dc').rglob('*') if p.is_file()])


def jvp_one(job):
    cfg=a.read(Path(job['config']))
    dest=LOCAL/'jvp'/job['kind']/job['model']/job['case']/f"vg_{job['index']:03d}"
    cfg.update(simulation_type='newton_jvp_probe',state_file=job['state'],output_csv=str(dest/'jvp.csv'))
    cfg.pop('output_state_file',None)
    if job['kind']=='post_solve':
        geo,mask=V.m.previous.prior.support(job)
        cfg['directions']=[dict(name=f'{mode}_{h:.0e}',mode=mode,amplitude_V=h,node_ids=np.where(mask)[0][::3].tolist(),exclude_contacts=True)
                           for mode in ('psi','phin','phip') for h in STEPS]
    v.write_same(dest/'config.json',cfg)
    if not (dest/'freeze.json').exists():d.matrix.freeze(dest/'freeze.json',[dest/'config.json',Path(job['state']),OUT/'validation_freeze.json'])
    a.verify(dest/'freeze.json')
    status=V.execute(dest/'config.json',RUNNER,V.environment());assert status['exit_code']==0,status
    rows=[]
    for r in a.rows(dest/'jvp.csv'):
        for block in ('psi','phin','phip'):
            ana,fd=float(r[f'analytic_{block}_norm']),float(r[f'finite_difference_{block}_norm'])
            absolute=float(r[f'{block}_relative_error'])*max(1.,fd)
            relative=absolute/max(ana,fd,1e-300)
            weak=(r['mode'],block) in (('phin','phip'),('phip','phin'))
            rows.append(dict(case=job['case'],device=job['device'],vd=job['vd'],index=job['index'],kind=job['kind'],model=job['model'],
                input_block=r['mode'],output_block=block,amplitude_V=float(r['amplitude_V']),analytic_norm=ana,fd_norm=fd,
                true_relative=relative,absolute_error=absolute,weak_carrier_cross_block=weak,qualified=relative<=1e-4))
    return rows


def jvp():
    a.verify(OUT/'dc_evidence.json')
    rows=a.rows(OUT/'attempts.csv');jobs=[];missing=[]
    for c in a.read(b.OUT/'contract.json')['cases']:
        for model,indices in (('old_slotboom',w.n.INDICES),('no_bgn',(40,50))):
            for index in indices:
                choices=[r for r in rows if r['case']==c['case'] and r['model']==model and int(r['index'])==index and r['qualified']=='True']
                if not choices:
                    missing.append(dict(case=c['case'],model=model,index=index));continue
                src=Path(next((r for r in choices if r['arm']=='vela'),choices[0])['dest'])
                jobs.append(dict(case=c['case'],device=c['device'],vd=c['vd'],model=model,index=index,kind='post_solve',config=str(src/'config.json'),state=str(src/'state.csv')))
        cfg=old_jvp.LOCAL/c['case']/'config.json'
        jobs.append(dict(case=c['case'],device=c['device'],vd=c['vd'],model='old_slotboom',index=40,kind='same_state',config=str(cfg),state=a.read(cfg)['state_file']))
    a.write(OUT/'jvp_execution.json',dict(jobs=jobs,missing=missing))
    with ThreadPoolExecutor(max_workers=2) as pool:result=[r for block in pool.map(jvp_one,jobs) for r in block]
    v.csv_union(OUT/'jvp_blocks.csv',result)
    d.matrix.freeze(OUT/'jvp_evidence.json',[OUT/'dc_evidence.json',OUT/'jvp_execution.json',OUT/'jvp_blocks.csv',old_jvp.OUT/'evidence.json']+[p for p in (LOCAL/'jvp').rglob('*') if p.is_file()])


def analyze():
    a.verify(OUT/'jvp_evidence.json')
    attempts=a.rows(OUT/'attempts.csv');native=a.rows(w.OUT/'native_points.csv')
    rows=[]
    for c in a.read(b.OUT/'contract.json')['cases']:
        geo,mask=V.m.previous.prior.support(c)
        for model,indices in (('old_slotboom',w.n.INDICES),('no_bgn',(40,50))):
            for index in indices:
                ref=next(r for r in native if r['case']==c['case'] and r['model']==model and int(r['index'])==index)
                row=dict(case=c['case'],device=c['device'],vd=c['vd'],vg=w.n.prior.GRID[index],index=index,model=model,
                         sentaurus_Id_A_per_um=float(ref['Id_A_per_um']),sentaurus_qualified=ref['native_qualified']=='True')
                chosen={}
                for arm in ('vela','native'):
                    group=[r for r in attempts if r['case']==c['case'] and r['model']==model and int(r['index'])==index and r['arm']==arm]
                    chosen[arm]=next((r for r in group if r['qualified']=='True'),group[-1])
                    row[arm+'_qualified']=chosen[arm]['qualified']=='True'
                    row[arm+'_Id_A_per_um']=float(chosen[arm].get('current_A_per_um') or 'nan')
                    row[arm+'_error_percent']=100*(row[arm+'_Id_A_per_um']/row['sentaurus_Id_A_per_um']-1)
                delta=dict(psi_max_V=math.inf,phin_max_V=math.inf,phip_max_V=math.inf,density_max_relative=math.inf)
                if all((Path(r['dest'])/'state.csv').exists() for r in chosen.values()):
                    delta=V.m.previous.prior.delta_states(*[d.ordered(Path(r['dest'])/'state.csv',geo.count) for r in chosen.values()],mask)
                error=abs(row['vela_Id_A_per_um']/row['native_Id_A_per_um']-1)
                row.update(**delta,dual_Id_relative=error)
                row['comparison_qualified']=v.dual_qualified(row['vela_qualified'],row['native_qualified'],delta,error) and row['sentaurus_qualified']
                baseline=a.rows((b.OUT if model=='old_slotboom' else w.OUT)/'comparison.csv')
                before=next(r for r in baseline if r['case']==c['case'] and int(r['index'])==index and r.get('model',model)==model)
                row['Id_change_from_pre_repair_relative']=row['vela_Id_A_per_um']/float(before['vela_Id_A_per_um'])-1
                rows.append(row)
    v.csv_union(OUT/'comparison.csv',rows)
    checks=a.rows(OUT/'jvp_blocks.csv')
    gated=[r for r in checks if r['weak_carrier_cross_block']=='False' and float(r['amplitude_V'])<=1e-5]
    summary=dict(points=len(rows),comparison_qualified=sum(r['comparison_qualified'] for r in rows),attempts=len(attempts),
        failed_attempts=sum(r['qualified']!='True' for r in attempts),jvp_block_checks=len(checks),jvp_gated_checks=len(gated),
        jvp_gated_failures=sum(r['qualified']!='True' for r in gated),jvp_gated_max_relative=max(float(r['true_relative']) for r in gated),
        max_Id_change_from_pre_repair_relative=max(abs(r['Id_change_from_pre_repair_relative']) for r in rows),
        numerical_repair_gate_passed=all(r['comparison_qualified'] for r in rows) and all(r['qualified']=='True' for r in gated))
    a.write(OUT/'summary.json',summary)
    d.matrix.freeze(OUT/'comparison_evidence.json',[OUT/'jvp_evidence.json',OUT/'comparison.csv',OUT/'summary.json'])
    print(summary,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('prepare','solve','jvp','analyze'))
    globals()[p.parse_args().action]()
