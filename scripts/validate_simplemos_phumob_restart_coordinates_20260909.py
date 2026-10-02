"""Isolate cancellation-free warm-restart coordinates on current production."""
import argparse
import copy
import shlex
import subprocess
from pathlib import Path
import calibrate_simplemos_phumob_local_source_20260909 as s

LOCAL=s.LOCAL/'restart_coordinates';OUT=s.OUT/'restart_coordinates';RUNNER=LOCAL/'runner.exe'

def build():
    LOCAL.mkdir(parents=True,exist_ok=False)
    source=s.REPO/'src/solver/NewtonSolver.cpp';text=source.read_text(encoding='utf-8')
    for car,short in (('electron','phin'),('hole','phip')):
        old=f'''static_cast<long double>(
                        initial.{car}QuasiFermiReferenceAt(i)) +
                    static_cast<long double>(initial.{short}Increment(i)) -
                    static_cast<long double>(
                        assembler.{car}QuasiFermiReferenceAt(node))'''
        new=f'''(static_cast<long double>(
                        initial.{car}QuasiFermiReferenceAt(i)) -
                     static_cast<long double>(
                        assembler.{car}QuasiFermiReferenceAt(node))) +
                    static_cast<long double>(initial.{short}Increment(i))'''
        assert text.count(old)==1;text=text.replace(old,new)
    target=LOCAL/'NewtonSolver.cpp';target.write_text(text,encoding='utf-8')
    entry=next(x for x in s.a.read(s.REPO/'build-release/compile_commands.json') if x['file'].endswith('/NewtonSolver.cpp'))
    cmd=shlex.split(entry['command'].replace('\\','/'));cmd=[('-DVELA_VERSION="0.1.0"' if x.startswith('-DVELA_VERSION=') else x) for x in cmd]
    cmd[cmd.index('-c')+1]=str(target);cmd[cmd.index('-o')+1]=str(LOCAL/'newton.o');s.a.write(LOCAL/'compile_command.json',cmd)
    r=subprocess.run(cmd,cwd=entry['directory'],env=s.V.environment(),capture_output=True,text=True);(LOCAL/'compile.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr[-3000:]
    link=s.a.read(s.LOCAL/'link_command.json');link[1]=str(LOCAL/'newton.o');link[link.index('-o')+1]=str(RUNNER);s.a.write(LOCAL/'link_command.json',link)
    r=subprocess.run(link,env=s.V.environment(),capture_output=True,text=True);(LOCAL/'link.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr[-3000:]
    s.d.matrix.freeze(OUT/'build_evidence.json',[Path(__file__).resolve(),source,s.prior.q.run.RUNNER,s.REPO/'build-release/libvela_core.a',s.OUT/'failed_trace/evidence.json']+[p for p in LOCAL.glob('*') if p.is_file()])

def prepare():
    s.a.verify(OUT/'build_evidence.json');jobs=[];files=[OUT/'build_evidence.json']
    case='m65_n23_vd_1p000000_endpoint';root=s.prior.LOCAL/'candidate/dc/phumob'/case
    bases=[('native_first',root/'native/vg_010/attempt_0'),('native_reload',root/'native/vg_010/attempt_1'),('qualified_control',root/'vela/vg_010/attempt_0')]
    for c in s.a.read(s.OUT/'contract.json')['cases']:bases.append((c['key'],Path(c['baseline'])))
    for label,base in bases:
        cfg=s.a.read(base/'config.json');device='n23' if 'n23' in str(base) else 'n19'
        for mode in ('legacy','corrected'):
            dest=LOCAL/'identity'/label/mode;deck=copy.deepcopy(cfg)
            deck.update(state_file=str(base/'state.csv'),output_state_file=str(dest/'state.csv'))
            deck['solver']['max_iter']=0;deck['solver']['carrier_row_convergence']['min_newton_max_iter']=0
            s.a.write(dest/'config.json',deck);s.V.post_config(deck,dest);jobs.append(dict(label=label,device=device,base=str(base),mode=mode,dest=str(dest)));files += [base/'state.csv',base/'all_row.csv']+list(dest.glob('*.json'))
    s.a.write(OUT/'contract.json',dict(jobs=jobs,scope='Zero-iteration pack/restart identity on 2 failed states, qualified same-bias control, and four Vg=0 controls. No source perturbation. Followed by at most one 200-iteration repair run from each failed output; keep all original failures.',gates=dict(state_relative_increment=1e-14,residual_identity_relative=1e-10,row=1e-6,kcl=1e-8),change='Subtract old-new references before adding saved increment in the two existing warm_start branches. Contact projection and consistency guards retained.',production_modified=False))
    s.d.matrix.freeze(OUT/'freeze.json',files+[OUT/'contract.json'])

def identity():
    s.a.verify(OUT/'freeze.json');rows=[]
    for j in s.a.read(OUT/'contract.json')['jobs']:
        dest=Path(j['dest']);runner=RUNNER if j['mode']=='corrected' else s.prior.q.run.RUNNER
        st=s.V.execute(dest/'config.json',runner,s.V.environment());assert st.get('iterations')==0,st
        s.V.execute(dest/'all_row.json',runner,s.V.environment())
        geo,mask=s.fields.support(j);old=s.d.ordered(Path(j['base'])/'all_row.csv',geo.count);new=s.d.ordered(dest/'all_row.csv',geo.count)
        residual=max(abs(float(x[k])-float(y[k]))/max(abs(float(x[car+'_flux_abs_sum'])),abs(float(x[car+'_recombination'])),1e-300) for x,y,keep in zip(old,new,mask) if keep for car,k in (('electron','electron_residual'),('hole','hole_residual')))
        before=s.d.ordered(Path(j['base'])/'state.csv',geo.count);after=s.d.ordered(dest/'state.csv',geo.count)
        delta=s.fields.delta_states(after,before,mask)
        increment=max(abs(float(x[car+'_qf_increment_V'])-float(y[car+'_qf_increment_V']))/max(abs(float(y[car+'_qf_increment_V'])),1e-300) for x,y,keep in zip(after,before,mask) if keep for car in ('electron','hole'))
        rows.append(dict(label=j['label'],mode=j['mode'],max_residual_change_over_row_scale=residual,max_increment_relative_change=increment,**delta,qualified=residual<=1e-10 and increment<=1e-14));print(rows[-1],flush=True)
    s.a.write_csv(OUT/'identity.csv',rows);s.d.matrix.freeze(OUT/'identity_evidence.json',[OUT/'freeze.json',OUT/'identity.csv']+[p for p in (LOCAL/'identity').rglob('*') if p.is_file()])

def solve():
    s.a.verify(OUT/'identity_evidence.json');assert all(r['qualified']=='True' for r in s.a.rows(OUT/'identity.csv') if r['mode']=='corrected')
    rows=[]
    for job in s.a.read(OUT/'contract.json')['jobs']:
        if job['mode']!='corrected' or not job['label'].startswith('native_'):continue
        base=Path(job['base']);dest=LOCAL/'dc'/job['label'];cfg=s.a.read(base/'config.json');cfg.update(state_file=str(base/'state.csv'),output_state_file=str(dest/'state.csv'))
        cfg['solver']['local_update_diagnostics']=dict(enabled=True,nodes=[967,983,1089],csv_file=str(dest/'updates.csv'),first_iterations=200,every_iterations=1)
        s.a.write(dest/'config.json',cfg);s.V.post_config(cfg,dest);s.d.matrix.freeze(dest/'freeze.json',[base/'state.csv',dest/'config.json',OUT/'identity_evidence.json'])
        st=s.V.execute(dest/'config.json',RUNNER,s.V.environment())
        for n in ('all_row','acceptance_edges'):s.V.execute(dest/(n+'.json'),RUNNER,s.V.environment())
        row=dict(label=job['label'],elapsed_seconds=st['elapsed_seconds'],**s.prior.q.run.w.old.prior.old.qualify(job,dest));rows.append(row);print(row,flush=True)
    s.a.write_csv(OUT/'dc.csv',rows);s.d.matrix.freeze(OUT/'dc_evidence.json',[OUT/'identity_evidence.json',OUT/'dc.csv']+[p for p in (LOCAL/'dc').rglob('*') if p.is_file()])

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('build','prepare','identity','solve'));globals()[p.parse_args().action]()
