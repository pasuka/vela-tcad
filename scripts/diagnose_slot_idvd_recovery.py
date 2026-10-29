"""Serial, bounded Fermi IdVd recovery diagnostics; never modifies a campaign ledger.

The source is a completed failed attempt. Local states and Poisson intermediates
remain diagnostic. Only the unchanged complete Newton gates can qualify a state.
"""
from pathlib import Path
import argparse, copy, csv, json, math, os, shutil, subprocess, time
from run_slot_ldmos_idl import qualify, sha, save
from run_slot_ldmos_noavalanche_continuation import redirect_diagnostics
from slot_local_row_repair import row_table
import slot_state_archive as states

ROOT=Path(__file__).resolve().parents[1]

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--runner',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--budget',type=float,default=1800.)
    p.add_argument('--child-cap',type=float,default=180.)
    a=p.parse_args()
    if not all(math.isfinite(x) and x>0 for x in (a.budget,a.child_cap)):p.error('Positive budgets required')
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    source=a.source.resolve();base=json.loads((source/'config.json').read_text())
    if base['solver']['carrier_statistics']!='fermi_dirac':raise ValueError('Explicit Fermi configuration required')
    seed=list((source/'rejected_states').glob('*_final.h5'))
    if len(seed)!=1:raise ValueError('One complete rejected final state required')
    seed=seed[0];bias=base['sweep']['bias_points'][0]
    runner=out/'runner.exe';shutil.copy2(a.runner,runner)
    shutil.copy2(__file__,out/'diagnose.py')
    protected={str(f):sha(f) for f in [source/'config.json',source/'iv.csv',seed,a.runner.resolve()]}
    protected.update({str(ROOT/base[k]):sha(ROOT/base[k]) for k in ('mesh_file','materials_file','node_doping_file')})
    summary=dict(completed=False,qualified=False,cases=[],budget_s=a.budget,child_cap_s=a.child_cap,
                 source=str(source),seed_sha256=sha(seed),runner_sha256=sha(runner),protected=protected,
                 policy=dict(max_dofs=64,alphas=[1.,.5],max_local_updates_per_arm=2,ledger_updates=0))
    env=os.environ.copy();env['PATH']='D:/msys64/ucrt64/bin;D:/msys64/usr/bin;'+env['PATH']
    started=time.monotonic()
    def persist():
        summary['elapsed_s']=time.monotonic()-started;save(out/'summary.json',summary)
    def execute(name,state,kind,**options):
        left=a.budget-(time.monotonic()-started)
        if left<2:raise TimeoutError('batch_budget_exhausted')
        d=out/name;d.mkdir();c=copy.deepcopy(base);redirect_diagnostics(c,d)
        c['solver']['performance_profiling']=dict(enabled=True,json_file=str(d/'performance.json'))
        c['output_csv']=str(d/'iv.csv')
        if kind=='strict' or kind=='poisson':
            c['sweep'].update(initial_state_file=str(state),write_state_file=str(d/'state.h5'))
            if kind=='strict':
                c['solver'].update(progress_directory=str(d/'progress'),soft_time_limit_seconds=.8*min(left-1,a.child_cap))
            else:
                c['solver'].update(method='poisson_only',poisson_only_quasi_fermi_initialization='preserve_initial')
        else:
            c.pop('sweep');c.update(simulation_type=kind,state_file=str(state))
            for contact in c['contacts']:
                if contact['name']=='drain':contact['bias']=bias
        c.update(options);save(d/'config.json',c)
        tick=time.monotonic();cap=min(a.child_cap,left-1)
        summary['active']=name;persist()
        with (d/'stdout.log').open('wb') as stdout,(d/'stderr.log').open('wb') as stderr:
            try:rc=subprocess.run([str(runner),'--config',str(d/'config.json')],cwd=ROOT,env=env,
                                   stdout=stdout,stderr=stderr,timeout=cap).returncode
            except subprocess.TimeoutExpired:rc=124
        entry=dict(name=name,returncode=rc,elapsed_s=time.monotonic()-tick,config_sha256=sha(d/'config.json'),qualified=False)
        if kind in ('strict','poisson') and (d/'iv.csv').exists():
            with (d/'iv.csv').open() as f:r=list(csv.DictReader(f))
            if r:entry['iv']=r[-1]
            entry['qualified']=kind=='strict' and rc==0 and (d/'state.h5').exists() and qualify(entry.get('iv',{}),bias)
            if entry['qualified']:entry['state_sha256']=sha(d/'state.h5')
        summary['cases'].append(entry);summary.pop('active',None);save(d/'exit.json',entry);persist()
        print(json.dumps(entry),flush=True)
        return d,entry
    def audit(name,state):
        d,e=execute(name+'_rows',state,'newton_carrier_term_probe',
                    carrier_term_probe=dict(solved_equation_terms=True),output_csv=str(out/(name+'_rows')/'terms.csv'))
        if e['returncode']:raise RuntimeError('row_probe_failed')
        rows=row_table(d/'terms.csv',base['solver']['carrier_row_convergence']['scale_floor'])
        g,ge=execute(name+'_global',state,'newton_candidate_conservation_probe',candidate_port_diagnostics=True,
                     output_csv=str(out/(name+'_global')/'result.json'))
        if ge['returncode']:raise RuntimeError('global_probe_failed')
        result=json.loads((g/'result.json').read_text())
        value=dict(max_ratio=max(r['ratio'] for r in rows.values()),bad_rows=sum(r['ratio']>.001 for r in rows.values()),
                   active_rows=len(rows),global_result=result)
        save(out/(name+'_audit.json'),value)
        return rows,value
    try:
        rows,baseline=audit('baseline',seed)
        mesh=json.loads(Path(base['mesh_file']).read_text());neighbors={}
        for tri in mesh['triangles']:
            ids=tri['node_ids'] if isinstance(tri,dict) else tri
            for i in ids:neighbors.setdefault(i,set()).update(j for j in ids if j!=i)
        bad=sorted((k for k,r in rows.items() if r['ratio']>.001),key=lambda k:-rows[k]['ratio'])
        if len(bad)>64:raise RuntimeError('patch_budget_exceeded')
        edge={(j,c) for i,c in bad for j in neighbors[i] if (j,c) in rows}-set(bad)
        patch=bad+sorted(edge,key=lambda k:-rows[k]['ratio'])[:64-len(bad)]
        save(out/'patch.json',dict(bad=bad,patch=patch,neighbor_rule='highest current ratio among active one-ring boundary'))
        summary['arms']=[]
        for alpha in (1.,.5):
            current=seed;arm=dict(alpha=alpha,steps=[]);summary['arms'].append(arm)
            for cycle in range(2):
                name=f'local_{alpha}_{cycle}'
                d,e=execute(name,current,'newton_block_step_probe',block_modes=['carrier_local'],
                      block_step_damping=alpha,local_carrier_nodes=[dict(node_id=i,carrier=c) for i,c in patch],
                      output_state_prefix=str(out/name/'trial'))
                if e['returncode']:arm['stop_reason']='local_probe_failed';break
                trial=d/'trial_carrier_local.h5';rr,av=audit(name,trial)
                arm['steps'].append(dict(state=str(trial),**av))
                # Diagnostic only; neither row descent nor this screen grants acceptance.
                if av['max_ratio']>= (baseline if cycle==0 else arm['steps'][-2])['max_ratio']:
                    arm['stop_reason']='no_all_row_descent';break
                current=trial
            arm['last_state']=str(current);persist()
        best=min(summary['arms'],key=lambda arm:arm['steps'][-1]['max_ratio'] if arm['steps'] else float('inf'))
        state=Path(best['last_state'])
        pd,pe=execute('poisson_repair',state,'poisson')
        if pe['returncode']==0 and (pd/'state.h5').exists() and pe.get('iv',{}).get('converged')=='1':state=pd/'state.h5'
        _,strict=execute('strict_after_repair',state,'strict');summary['qualified']=strict['qualified']
        # Separate control: no local or Poisson repair, same rejected seed.
        execute('strict_control',seed,'strict')
        # Three physical amplitudes plus 1 ULP in physical QF. A null FD response
        # is evidence about representation, never evidence of exact derivatives.
        fields,_=states.read(seed,base,ROOT)
        focus=[i for i,c in bad if c=='hole'][:3]
        if not focus:raise RuntimeError('no_weak_hole_rows_for_targeted_diagnostic')
        for index,h in enumerate([1e-5,1e-6,1e-7,math.ulp(float(fields['phip'][focus[0]]))]):
            name=f'jvp_{index}'
            execute(name,seed,'newton_jvp_probe',directions=[dict(name='weak_holes',mode='phip',amplitude_V=h,
                    node_ids=focus,node_index_base=0,exclude_contacts=True)],
                    sample_rows=[dict(node_id=i,block='phip') for i in focus],row_output_csv=str(out/name/'samples.csv'))
        summary['completed']=True
    except (TimeoutError,RuntimeError) as exc:
        summary.update(completed=True,stop_reason=str(exc))
    finally:
        summary['protected_unchanged']=all(sha(Path(f))==h for f,h in protected.items())
        persist()
    return 0 if summary.get('completed') and summary['protected_unchanged'] else 2

if __name__=='__main__':raise SystemExit(main())
