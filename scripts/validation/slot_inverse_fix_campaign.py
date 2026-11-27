"""Finite same-parent mapping trials, repeat checks and prospective P2 samples.

Only invokes isolated runners. No long sweep or acceptance-rule changes.
"""
from pathlib import Path
import copy, datetime, json, math, os, shutil, sys, time, traceback
import h5py, numpy as np

R=Path.cwd(); B=R.parent
sys.path.insert(0,str(B/'s21p1_parallel'))
import common as C
C.END=datetime.datetime.fromisoformat('2026-10-22T19:16:00+08:00').timestamp()
OLD=B/'s22p2_scans12_1103/idvd'; FAIL=OLD/'v2_point_0035_0.155955974V'
NEW=R/'build/runner'; BASE=C.get(FAIL/'out/config.json')
E={k:v for k,v in os.environ.items() if not k.startswith('VELA_')}
E.update(C.get(OLD/'environment.json'))
S=dict(completed=False,phase='preflight',deadline_epoch=C.END,records=[],comparisons=[],
       acceptance_rules_changed=False,p2_thresholds='not_calibrated')
def save(): C.put(R/'status.json',S)
def config(d,parent,target,enabled,updates=25):
    c=copy.deepcopy(BASE); c['state_file']=str(parent)
    c['output_state_file']=str(d/'state.h5');c['output_csv']=str(d/'rows.csv')
    for contact in c['contacts']:
        if contact['name']=='drain':contact['bias']=target
    c['solver']['progress_directory']=str(d/'progress')
    c['solver']['carrier_row_convergence']['diagnostic_csv']=str(d/'qualification.csv')
    c['population_memory'].update(output_directory=str(d/'out'),max_updates=updates,
        soft_seconds=900,validation_snapshots=True,validation_snapshot_max_updates=25,
        capture_steps=[],update_diagnostics=True,projection_diagnostics=True,
        dilute_log_inverse=enabled)
    return c
def environment(d):
    e=E.copy();e.update({k:'1' for k in ['VELA_P1_COMPLETE_CAPTURE','VELA_P1_SLOT_STAGE_CAPTURE',
        'VELA_P1_SLOT_GUARD_CAPTURE','VELA_P1_GUARD_CAPTURE','VELA_P1_SHADOW_CAPTURE']})
    e.update(VELA_P1_ENTRY_CAPTURE_DIR=str(d/'entries'),
        VELA_P1_CONSTITUTIVE_DIR=str(d/'constitutive'),
        VELA_P1_CONSTITUTIVE_MAX_CALLS='32',VELA_P1_CONSTITUTIVE_STATE_TAG=d.name)
    return e
def parsed(p):
    t=Path(p).read_text();return json.JSONDecoder().raw_decode(t[t.find('{'):])[0]
def promotable(q,rc):
    return (rc==0 and q.get('qualified') is True and q.get('converged') is True
        and q.get('strict_qualification_satisfied') is True
        and q.get('stop_reason')=='SLOT_rhs_or_raw_update'
        and q.get('readback_audit',{}).get('consistent') is True
        and (q.get('rhs_satisfied') is True or q.get('raw_update_satisfied') is True))
def run(name,parent,target,enabled,updates=25,binary=NEW,prepare=False,mechanism=False):
    d=R/name;d.mkdir();S['phase']=name;save()
    rec=dict(name=name,parent=str(parent),parent_sha256=C.sha(parent),target_V=target,
             enabled=enabled,mechanism_only=mechanism,promoted=False)
    input_parent=parent
    try:
        if prepare:
            p=d/'prepare';p.mkdir();c=config(p,parent,target,enabled)
            c.pop('population_memory');c['simulation_type']='newton_zero_update_qualification_probe'
            c['solver']['max_iter']=0;c['solver']['global_continuity_closure']['feasibility_restoration']='off'
            C.put(p/'config.json',c)
            ex=C.run(p,'prepare',[binary,'--config',p/'config.json'],env=E,cap=300,check=False)
            q=parsed(p/'prepare.stdout');C.put(p/'result.json',q)
            assert ex['exit'] in (0,1) and q.get('iterations')==0 and math.isfinite(q['final_residual']) and (p/'state.h5').exists()
            parent=p/'state.h5'
            b=d/'boundary';b.mkdir();c=config(b,parent,target,enabled)
            c.pop('population_memory');c['simulation_type']='newton_poisson_term_probe'
            C.put(b/'config.json',c);C.run(b,'boundary',[binary,'--config',b/'config.json'],env=E,cap=240)
            import csv
            rows=list(csv.DictReader((b/'rows.csv').open()))
            cons=[x for x in rows if x['constrained'] in ('1','true')]
            error=max((abs(float(x['production_residual'])) for x in cons),default=float('inf'))
            assert error<1e-10,'target boundary preparation failed'
            rec['initial_boundary_max']=error
        c=config(d,parent,target,enabled,updates);C.put(d/'config.json',c)
        env=environment(d);C.put(d/'environment.json',env)
        ex=C.run(d,'solve',[binary,'--config',d/'config.json'],env=env,cap=1050,check=False)
        q=C.get(d/'out/result.json') if (d/'out/result.json').exists() else {}
        hist=q.get('history',[])
        rec.update(exit=ex['exit'],stop_reason=q.get('stop_reason'),error=q.get('error'),
            updates=len(hist),converged=q.get('converged'),qualified=q.get('qualified'),
            final_residual=q.get('final_residual'),current=q.get('contact_currents_A_per_um'),
            readback=q.get('readback_audit'),wall_s=ex['wall_s'],
            log_selected=sum(x.get('finite_update_arithmetic',{}).get('dilute_log_selected',0) for x in hist),
            log_applied=sum(x.get('finite_update_arithmetic',{}).get('dilute_log_applied',0) for x in hist))
        if (d/'out/inverse_failure.json').exists():rec['inverse_failure']=C.get(d/'out/inverse_failure.json')
        rec['completed_attempt']=ex['exit']!=124
        state=d/'out/final.h5'
        if not mechanism and promotable(q,ex['exit']) and state.exists():
            with h5py.File(state) as f:
                assert int(f.attrs['node_count'])==13881
                md=json.loads(f.attrs['metadata_json']);assert abs(md['contact_biases_V']['drain']-target)<1e-14
                assert all(np.isfinite(f['fields'][key][:]).all() for key in f['fields'])
            rec.update(promoted=True,state=str(state),state_sha256=C.sha(state))
        assert C.sha(input_parent)==rec['parent_sha256']
    except TimeoutError:
        rec.update(error='deadline',completed_attempt=False);raise
    except Exception as e:rec.update(error=repr(e),traceback=traceback.format_exc())
    finally:
        C.put(d/'record.json',rec);S['records'].append(rec);save()
    return rec
def compare(a,b):
    qa=C.get(R/a/'out/result.json');qb=C.get(R/b/'out/result.json')
    z=dict(a=a,b=b,history_steps=[len(qa.get('history',[])),len(qb.get('history',[]))])
    keys=['rho','raw_error','projected','psi_damped','rows','closure','boundary_preserved',
          'mapping_error','population_extrema','contact_currents_A_per_um']
    z['core_history_equal']=([{k:x.get(k) for k in keys} for x in qa.get('history',[])]==
                             [{k:x.get(k) for k in keys} for x in qb.get('history',[])])
    z['currents_equal']=qa.get('contact_currents_A_per_um')==qb.get('contact_currents_A_per_um')
    for name in ['final.h5','failed.h5']:
        pa=R/a/'out'/name;pb=R/b/'out'/name
        if pa.exists() and pb.exists():
            with h5py.File(pa) as f,h5py.File(pb) as g:
                z['numeric_fields_equal']=set(f['fields'])==set(g['fields']) and all(
                    f['fields'][k][:].tobytes()==g['fields'][k][:].tobytes() for k in f['fields'])
            break
    S['comparisons'].append(z);save();return z
def main():
    (R/'execution.lock').mkdir()
    assert C.get(R/'build_status.json')['passed'] is True
    manifest=C.get(R/'manifest.json')
    assert all(C.sha(p)==h for p,h in manifest.items())
    inputs=[Path(BASE[k]) for k in ['mesh_file','materials_file','node_doping_file','state_file']]
    inputs += [Path(BASE['population_memory']['support_file']), FAIL/'out/failed_matrix_parent.h5',OLD/'inputs/parent6mV.h5']
    bound={str(p):C.sha(p) for p in inputs}
    C.put(R/'contract.json',dict(deadline_epoch=C.END,input_hashes=bound,backend=E['VELA_LINEAR_SOLVER'],
        threads=1,gate_V=5.,mesh_nodes=13881,source_campaign=str(OLD),
        cases=['frozen-step off/on','full-target repeats','prospective7.5mV','prospective8mV'],
        prospective_samples='Selected before this build executes; mechanism samples not held-out',
        acceptance='legacy unchanged; incomplete and diagnostic states never promoted',
        p2='shadow only; no threshold calibration; missing fields not_evaluated'))
    parent=FAIL/'out/failed_matrix_parent.h5';target=.15595597396240235
    run('mechanism_old',parent,target,False,1,OLD/'bin/runner',mechanism=True)
    run('mechanism_off',parent,target,False,1,mechanism=True)
    cmp=compare('mechanism_old','mechanism_off')
    assert cmp.get('numeric_fields_equal') and cmp['core_history_equal'] and cmp['currents_equal'],'default-off invariance failed'
    for i in range(1,4):run('mechanism_on'+str(i),parent,target,True,1,mechanism=True)
    compare('mechanism_on1','mechanism_on2');compare('mechanism_on1','mechanism_on3')
    # Full target begins at the original prepared initial guess, not a failed candidate.
    run('full_off',Path(BASE['state_file']),target,False)
    for i in range(1,3):run('full_on'+str(i),Path(BASE['state_file']),target,True)
    compare('full_on1','full_on2')
    parent=OLD/'inputs/parent6mV.h5'
    run('prospective_7p5_off',parent,.0075,False,prepare=True)
    q=run('prospective_7p5_on',parent,.0075,True,prepare=True)
    if q.get('promoted'):
        run('adjacent_8_on',Path(q['state']),.008,True,prepare=True)
        run('direct_8_on',parent,.008,True,prepare=True)
        compare('adjacent_8_on','direct_8_on')
    S['hashes_unchanged']=all(C.sha(p)==h for p,h in {**bound,**manifest}.items())
    assert S['hashes_unchanged']
    S['phase']='completed'
if __name__=='__main__':
    try:main()
    except Exception as e:S.update(phase='deadline' if isinstance(e,TimeoutError) else 'controller_error',error=repr(e),traceback=traceback.format_exc())
    finally:S.update(completed=True,completed_epoch=time.time());save()
    sys.exit(1 if S.get('error') else 0)
