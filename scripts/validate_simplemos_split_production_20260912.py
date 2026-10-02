"""Opt-in split-state production migration; unchanged original 16-point cohort.
No native jobs or model restoration are launched by this driver.
"""
import argparse, copy, json
from pathlib import Path
import validate_simplemos_production_weighted_20260910 as previous

a,d,run=previous.a,previous.d,previous.run
REPO=previous.REPO
LOCAL=REPO/'build-release/split_prod_0912/v3'
OUT=REPO/'reference_tcad/simplemos_sentaurus2022/phumob_split_production_20260912/v3'
BASE=REPO/'reference_tcad/simplemos_sentaurus2022/phumob_packed_restart_20260910/cohort_v2/validation_contract.json'
SOURCES=('src/equation/SplitDDRuntime.cpp','include/vela/equation/SplitDDRuntime.h',
 'include/vela/equation/SplitDDOperator.h','include/vela/numerics/SplitDDState.h',
 'src/equation/CoupledDDAssembler.cpp','include/vela/equation/CoupledDDAssembler.h',
 'src/solver/NewtonSolver.cpp','include/vela/solver/NewtonSolver.h','include/vela/solver/GummelSolver.h',
 'src/numerics/LineSearch.cpp','include/vela/numerics/LineSearch.h','src/io/DDSolutionCsv.cpp',
 'src/post/ContactCurrent.cpp','include/vela/post/ContactCurrent.h','src/tools/vela_example_runner.cpp',
 'tests/test_split_dd_runtime.cpp','CMakeLists.txt')

def configure():
    previous.LOCAL=LOCAL/'cohort';previous.OUT=OUT;previous.OPTIONS=dict(split_dd_state=True)
    previous.configure()

def prepare():
    contract=copy.deepcopy(a.read(BASE)); files=[BASE,Path(__file__).resolve(),previous.RUNNER,REPO/'build-release/libvela_core.a']
    for job in contract['jobs']:
        old=Path(job['original_config']);cfg=a.read(old)
        assert cfg['solver']['mobility']['model']=='phumob' and cfg['solver']['mobility']['edge_averaging']=='element_box_phumob'
        cfg['solver']['split_dd_state']=True
        dest=LOCAL/'inputs'/job['case']/str(job['index'])/job['arm']/'config.json'
        a.write(dest,cfg);job['original_config']=str(dest)
        files += [old,dest,Path(job['seed'])]+[Path(cfg[k]) for k in ('mesh_file','materials_file','node_doping_file')]
    contract.update(change='Only opt in to the unified split-coordinate state evaluator and checkpoint. Preserve original seeds, physical parameters, encoded boundary/gauge equations, solver gates, step caps and one-reload limit.',options=dict(split_dd_state=True),
      staged_release=['fixed-state and derivative/restart/port preflight','16/16 dual initialization','four full PhuMob curves','Enormal','high-field saturation','original full cohort'],
      precision='100 decimal digits for state-dependent residual, SRH, PhuMob and port evaluation; existing analytic double Jacobian and SparseLU with four refinements remain subject to derivative qualification.',
      diagnostic_reference='Split-state *LongDoubleReference columns project the same wide evaluator; these columns are not independent accuracy evidence.')
    assert len(contract['jobs'])==32
    a.write(OUT/'validation_contract.json',contract)
    d.matrix.freeze(OUT/'validation_freeze.json',files+[REPO/p for p in SOURCES]+[OUT/'validation_contract.json'])
    print('Frozen unchanged 16 points, 32 original seeds and opt-in implementation.',flush=True)

def preflight():
    configure();a.verify(OUT/'validation_freeze.json');results=[];files=[OUT/'validation_freeze.json']
    case='m65_n23_vd_1p000000_endpoint'
    for arm in ('vela','native'):
        source=REPO/'build-release/pr_fix/cohort_v2/dc/phumob'/case/arm/'vg_000/attempt_0'
        cfg=a.read(source/'config.json');cfg['solver'].update(split_dd_state=True,max_iter=0)
        cfg['solver']['carrier_row_convergence']['mode']='report'
        dest=LOCAL/'preflight'/arm
        cfg.update(state_file=str(source/'state.csv'),output_state_file=str(dest/'state.csv'))
        a.write(dest/'config.json',cfg);first=run.V.execute(dest/'config.json',previous.RUNNER,run.V.environment())
        assert (dest/'state.csv').exists(),first
        reload=copy.deepcopy(cfg);reload.update(state_file=str(dest/'state.csv'),output_state_file=str(dest/'reload_state.csv'))
        a.write(dest/'reload.json',reload);again=run.V.execute(dest/'reload.json',previous.RUNNER,run.V.environment())
        identity=first['final_block_residuals']==again['final_block_residuals'] and (dest/'state.csv').read_bytes()==(dest/'reload_state.csv').read_bytes()
        run.V.post_config(cfg,dest)
        for name in ('all_row','acceptance_edges'):
            status=run.V.execute(dest/(name+'.json'),previous.RUNNER,run.V.environment());assert status['exit_code']==0,status
        terms=a.rows(dest/'all_row.csv');edges=a.rows(dest/'acceptance_edges.csv')
        import math
        from collections import defaultdict
        ef,hf=defaultdict(list),defaultdict(list)
        for e in edges:
            for key,store in (('electron_flux',ef),('hole_flux',hf)):
                store[int(e['node0'])].append(float(e[key]));store[int(e['node1'])].append(-float(e[key]))
        reconstruction=[]
        for r in terms:
            i=int(r['node_id'])
            for car,store in (('electron',ef),('hole',hf)):
                scale=float(r[car+'_flux_abs_sum'])
                if scale<=0:continue
                reconstruction.append(abs(math.fsum(store[i])-float(r[car+'_flux']))/max(scale,1e-300))
        assert reconstruction,'No active rows selected'
        fun=a.read(dest/'acceptance_edges.json');fun.pop('output_csv',None)
        fun.update(simulation_type='terminal_current_functional_probe',contact='drain',residual_output_csv=str(dest/'port_residual.csv'),contact_edge_output_csv=str(dest/'port_edges.csv'))
        a.write(dest/'functional.json',fun);f=run.V.execute(dest/'functional.json',previous.RUNNER,run.V.environment());assert f['exit_code']==0,f
        port_error=abs(f['current_A_per_um']/f['contact_current_extractor_A_per_um']-1)
        jvp=copy.deepcopy(fun);jvp.pop('residual_output_csv');jvp.pop('contact_edge_output_csv')
        jvp.update(simulation_type='newton_jvp_probe',output_csv=str(dest/'jvp.csv'),directions=[dict(name=f'{mode}_{h:.0e}',mode=mode,amplitude_V=h,node_ids=[1000,1009,1089],exclude_contacts=True) for mode in ('psi','phin','phip') for h in (1e-6,1e-20)])
        a.write(dest/'jvp.json',jvp);js=run.V.execute(dest/'jvp.json',previous.RUNNER,run.V.environment());assert js['exit_code']==0,js
        checks=[]
        for r in a.rows(dest/'jvp.csv'):
            for block in ('psi','phin','phip'):
                if (r['mode'],block) in (('phin','phip'),('phip','phin')):continue
                an=float(r[f'analytic_{block}_norm']);fd=float(r[f'finite_difference_{block}_norm'])
                checks.append(float(r[f'{block}_relative_error'])*max(1.,fd)/max(an,fd,1e-300))
        result=dict(arm=arm,restart_exact=identity,flux_reconstruction=max(reconstruction),port_relative=port_error,jvp_max=max(checks))
        result['qualified']=identity and max(reconstruction)<1e-12 and port_error<=1e-8 and max(checks)<=1e-4
        results.append(result);print(result,flush=True)
        files += [source/'config.json',source/'state.csv']+[p for p in dest.rglob('*') if p.is_file()]
    a.write(OUT/'preflight.json',results);d.matrix.freeze(OUT/'preflight_evidence.json',files+[OUT/'preflight.json'])
    assert all(r['qualified'] for r in results),results

def pilot():
    configure();a.verify(OUT/'validation_freeze.json')
    a.verify(OUT/'preflight_evidence.json');assert all(r['qualified'] for r in a.read(OUT/'preflight.json'))
    job=next(j for j in a.read(OUT/'validation_contract.json')['jobs'] if j['device']=='n23' and j['vd']==1 and j['index']==0 and j['arm']=='native')
    rows=run.solve_job(job)
    a.write(OUT/'pilot.json',rows)
    print(rows,flush=True)

def solve():
    configure();a.verify(OUT/'preflight_evidence.json');assert all(r['qualified'] for r in a.read(OUT/'preflight.json'))
    run.solve()
def analyze(): configure();previous.analyze()
def jvp(): configure();previous.jvp()

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('prepare','preflight','pilot','solve','analyze','jvp'));globals()[p.parse_args().action]()
