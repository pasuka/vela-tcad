"""Unified production numerics: unchanged sixteen PhuMob controls, dual seeds."""
import argparse
import copy
import math
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import validate_simplemos_phumob_lowvg_20260909 as previous

a,d,run=previous.a,previous.d,previous.q.run
REPO=previous.REPO
LOCAL=REPO/'build-release/phumob_numerics_production_20260910/weighted_merit'
OUT=REPO/'reference_tcad/simplemos_sentaurus2022/phumob_numerics_production_20260910/weighted_merit'
RUNNER=REPO/'build-release/vela_example_runner.exe'
OPTIONS=dict(poisson_residual_precision='binary128',stable_merit_comparison=True,exact_dirichlet_updates=True)
SOURCES=('src/solver/NewtonSolver.cpp','include/vela/solver/NewtonSolver.h',
         'src/equation/CoupledDDAssembler.cpp','include/vela/equation/CoupledDDAssembler.h',
         'include/vela/equation/ExtendedPoissonResidual.h','include/vela/numerics/StableMeritComparison.h')
NATIVES=(previous.NOUT/'native_points.csv',previous.comparison.v.previous.p.OUT/'supported_export/native_points.csv')


_post_config=run.V.post_config

def post_config(cfg,dest):
    probe=copy.deepcopy(cfg)
    probe['solver']['stable_merit_comparison']=False
    return _post_config(probe,dest)

def configure():
    run.LOCAL,run.OUT,run.RUNNER=LOCAL,OUT,RUNNER
    run.V.post_config=post_config


def prepare():
    a.verify(OUT.parent/'unit_block_guard/restricted_norm_evidence.json')
    paths=(previous.OUT/'contract.json',previous.previous.OUT/'contract.json')
    contracts=[a.read(p) for p in paths]
    jobs=[]; files=[*paths,*NATIVES,Path(__file__).resolve(),Path(run.__file__),RUNNER,REPO/'build-release/libvela_core.a',OUT.parent/'unit_block_guard/restricted_norm_evidence.json']
    files += [REPO/p for p in SOURCES]
    for old in [j for c in contracts for j in c['jobs']]:
        cfg=a.read(Path(old['original_config'])); new=copy.deepcopy(cfg)
        assert new['solver']['mobility']['edge_averaging']=='element_box_phumob'
        assert new['solver']['mobility']['model']=='phumob'
        assert new['solver']['global_continuity_closure']['mode']=='off'
        assert new['solver'].get('residual_norm','block') in ('l2','block')
        new['solver'].update(OPTIONS)
        dest=LOCAL/'inputs'/old['case']/f"vg_{old['index']:03d}"/old['arm']/'config.json'
        a.write(dest,new); jobs.append(dict(old,original_config=str(dest)))
        files += [dest,Path(old['original_config']),Path(old['seed'])]
        files += [Path(new[k]) for k in ('mesh_file','materials_file','node_doping_file')]
    assert len(jobs)==32 and len({(j['case'],j['index']) for j in jobs})==16
    contract=dict(jobs=jobs,points=16,initial_attempts=32,gates=contracts[0]['gates'],options=OPTIONS,
        physics_changed=False,acceptance_changed=False,independent_seeds_preserved=True,
        recovery='At most one same-bias reload following a failed attempt. Preserve every attempt.',
        jvp_gate=1e-4,weak_cross='Retained, not qualified from whole-residual differences.',
        source_relative_floor=1e-10,finite_SRH_volume_replacement=False)
    a.write(OUT/'validation_contract.json',contract)
    d.matrix.freeze(OUT/'validation_freeze.json',files+[OUT/'validation_contract.json'])
    print('Frozen production binary, 16 points and 32 original independent seeds.',flush=True)


def solve():
    configure();run.solve()


def analyze():
    configure();a.verify(OUT/'dc_evidence.json')
    attempts=a.rows(OUT/'attempts.csv'); native=[r for p in NATIVES for r in a.rows(p)]
    cases={(j['case'],j['index']):j for j in a.read(OUT/'validation_contract.json')['jobs']}
    points=[]; selected=[]; activation=[]
    for (case,index),job in cases.items():
        ref=next(r for r in native if r['case']==case and int(r['index'])==index and r['model']=='phumob')
        row=dict(case=case,device=job['device'],vd=job['vd'],vg=job['vg'],index=index,
                 sentaurus_Id_A_per_um=float(ref['Id_A_per_um']),sentaurus_qualified=ref['native_qualified']=='True')
        chosen=[]
        for arm in ('vela','native'):
            group=[r for r in attempts if r['case']==case and int(r['index'])==index and r['arm']==arm]
            r=next((r for r in group if r['qualified']=='True'),group[-1]); chosen.append(r);selected.append(r)
            current=float(r.get('current_A_per_um') or 'nan')
            row.update({arm+'_qualified':r['qualified']=='True',arm+'_Id_A_per_um':current,
                        arm+'_error_percent':100*(current/row['sentaurus_Id_A_per_um']-1)})
            p=Path(r['dest'])/'independent_acceptance.json'
            if p.exists():
                activation.extend(dict(case=case,index=index,arm=arm,carrier=k,**val) for k,val in a.read(p)['closure'].items())
        geo,mask=run.V.m.previous.prior.support(job)
        delta=dict(psi_max_V=math.inf,phin_max_V=math.inf,phip_max_V=math.inf,density_max_relative=math.inf)
        if all((Path(r['dest'])/'state.csv').exists() for r in chosen):
            delta=run.V.m.previous.prior.delta_states(*[d.ordered(Path(r['dest'])/'state.csv',geo.count) for r in chosen],mask)
        relative=abs(row['vela_Id_A_per_um']/row['native_Id_A_per_um']-1)
        row.update(**delta,dual_Id_relative=relative)
        row['comparison_qualified']=row['sentaurus_qualified'] and run.v.dual_qualified(row['vela_qualified'],row['native_qualified'],delta,relative)
        # Native solver qualification and the native-initialized Vela run are separate gates.
        row['comparison_qualified'] &= ref['native_qualified']=='True'
        points.append(row)
    run.v.csv_union(OUT/'comparison.csv',points);run.v.csv_union(OUT/'selected_states.csv',selected)
    run.v.csv_union(OUT/'source_activation.csv',activation)
    summary=dict(points=len(points),qualified=sum(r['comparison_qualified'] for r in points),attempts=len(attempts),
        failed_attempts=sum(r['qualified']!='True' for r in attempts),qualified_selected=sum(r['qualified']=='True' for r in selected),
        max_row_ratio=max(float(r.get('max_row_ratio') or 'inf') for r in selected),
        max_dual_Id_relative=max(r['dual_Id_relative'] for r in points),
        max_native_error_percent=max(abs(r['vela_error_percent']) for r in points),
        source_components=len(activation),source_relative_active=sum(r['qualified'] for r in activation),
        acceptance_changed=False,options=OPTIONS)
    a.write(OUT/'comparison_summary.json',summary)
    d.matrix.freeze(OUT/'comparison_evidence.json',[OUT/'dc_evidence.json',*NATIVES]+[OUT/n for n in ('comparison.csv','selected_states.csv','source_activation.csv','comparison_summary.json')])
    print(summary,flush=True)


def jvp():
    configure();a.verify(OUT/'comparison_evidence.json')
    jobs=[]
    for r in a.rows(OUT/'selected_states.csv'):
        if r['arm']!='vela' or r['qualified']!='True':continue
        p=Path(r['dest']);jobs.append(dict(case=r['case'],device=r['device'],vd=float(r['vd']),index=int(r['index']),
            model='phumob',kind='post_solve',config=str(p/'config.json'),state=str(p/'state.csv')))
    with ThreadPoolExecutor(max_workers=2) as pool:rows=[r for block in pool.map(run.jvp_one,jobs) for r in block]
    run.v.csv_union(OUT/'jvp_blocks.csv',rows)
    gated=[r for r in rows if not r['weak_carrier_cross_block'] and r['amplitude_V']<1e-4]
    summary=dict(states=len(jobs),checks=len(gated),failures=sum(not r['qualified'] for r in gated),max_relative=max(r['true_relative'] for r in gated))
    a.write(OUT/'jvp_summary.json',summary)
    d.matrix.freeze(OUT/'jvp_evidence.json',[OUT/'comparison_evidence.json',OUT/'jvp_blocks.csv',OUT/'jvp_summary.json']+[p for p in (LOCAL/'jvp').rglob('*') if p.is_file()])
    print(summary,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('prepare','solve','analyze','jvp'))
    globals()[p.parse_args().action]()
