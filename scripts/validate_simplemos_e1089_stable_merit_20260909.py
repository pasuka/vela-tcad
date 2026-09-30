"""Single-axis exact norm-difference isolation, with original convergence gates."""
import argparse
import copy
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import validate_simplemos_restart_source_followup_20260909 as w
import audit_simplemos_e1089_linesearch_20260909 as audit
import build_simplemos_stable_merit_20260908 as prior
s=w.s;LOCAL=w.LOCAL/'stable_merit';OUT=w.OUT/'stable_merit';RUNNER=LOCAL/'runner.exe'

def build():
    s.a.verify(audit.OUT/'precision/evidence.json');LOCAL.mkdir(parents=True,exist_ok=False)
    source=(s.LOCAL/'restart_coordinates/NewtonSolver.cpp').read_text(encoding='utf-8')
    marker='        const auto runLineSearch = [&](const VectorXd& trialStep) {';assert source.count(marker)==1;source=source.replace(marker,prior.PREFIX+'\n'+marker)
    start=source.index(marker);end=source.index('        auto ls = runLineSearch(step);',start);frag=source[start:end]
    target='                    : BacktrackingLineSearch::DecreaseAcceptFunction{});';assert frag.count(target)==1;frag=frag.replace(target,'                    : (stableMerit ? BacktrackingLineSearch::DecreaseAcceptFunction(stableDecrease) : BacktrackingLineSearch::DecreaseAcceptFunction{}));')
    (LOCAL/'NewtonSolver.cpp').write_text('#include "simplemos_stable_merit.hpp"\n#include <cstdlib>\n'+source[:start]+frag+source[end:],encoding='utf-8')
    cmd=s.a.read(s.LOCAL/'restart_coordinates/compile_command.json');cmd.insert(1,'-I'+str(prior.HEADER.parent));cmd[cmd.index('-c')+1]=str(LOCAL/'NewtonSolver.cpp');cmd[cmd.index('-o')+1]=str(LOCAL/'newton.o');s.a.write(LOCAL/'compile_command.json',cmd)
    r=subprocess.run(cmd,cwd=s.REPO/'build-release',env=s.V.environment(),capture_output=True,text=True);(LOCAL/'compile.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr[-3000:]
    link=s.a.read(w.LOCAL/'link_command.json');link[1]=str(LOCAL/'newton.o');link[link.index('-o')+1]=str(RUNNER);s.a.write(LOCAL/'link_command.json',link)
    r=subprocess.run(link,env=s.V.environment(),capture_output=True,text=True);(LOCAL/'link.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr[-3000:]
    test=(s.REPO/'tests/diagnostics/test_simplemos_stable_merit.cpp').read_text();assert test.count('count==26')==1 and test.count('descending==7')==1
    test=test.replace('26 frozen trial signs','13 electron-1089 trial signs').replace('count==26','count==13').replace('descending==7','descending==8');(LOCAL/'test.cpp').write_text(test)
    rows=s.a.rows(audit.OUT/'precision/merit.csv');root=audit.LOCAL/'failed_e1089/trace';fixtures=[]
    for i in range(13):
        row=next(r for r in rows if r['mode']=='production' and r['label']=='trial_'+str(i));fixtures.append(str(root/'base_direction.csv')+'\t'+str(root/('trial_'+str(i)+'_direction.csv'))+'\t'+('1' if row['would_decrease']=='True' else '0'))
    (LOCAL/'fixtures.tsv').write_text('\n'.join(fixtures)+'\n')
    testcmd=[cmd[0]]+[x for x in cmd[1:] if x.startswith(('-I','-D','-O','-std='))]+['-ID:/msys64/ucrt64/include/eigen3',str(LOCAL/'test.cpp'),'D:/msys64/ucrt64/lib/libCatch2Main.a','D:/msys64/ucrt64/lib/libCatch2.a','-o',str(LOCAL/'test.exe')]
    s.a.write(LOCAL/'test_compile.json',testcmd);r=subprocess.run(testcmd,env=s.V.environment(),capture_output=True,text=True);(LOCAL/'test_compile.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr[-3000:]
    env=s.V.environment();env['VELA_MERIT_FIXTURE_MANIFEST']=str(LOCAL/'fixtures.tsv');r=subprocess.run([str(LOCAL/'test.exe')],env=env,capture_output=True,text=True);(LOCAL/'test_result.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stdout+r.stderr
    seal();print(r.stdout,flush=True)

def seal():
    s.d.matrix.freeze(OUT/'build_evidence.json',[Path(__file__).resolve(),Path(prior.__file__),prior.HEADER,s.REPO/'tests/diagnostics/test_simplemos_stable_merit.cpp',audit.OUT/'precision/evidence.json',w.OUT/'freeze.json']+[p for p in LOCAL.iterdir() if p.is_file()])

def finish():
    # Preserve the first test command/log: its flag extraction omitted the
    # paired -isystem Eigen include. Solver compilation/linking had succeeded.
    cmd=s.a.read(LOCAL/'test_compile.json');cmd.insert(1,'-ID:/msys64/ucrt64/include/eigen3');s.a.write(LOCAL/'test_compile_resolved.json',cmd)
    r=subprocess.run(cmd,env=s.V.environment(),capture_output=True,text=True);(LOCAL/'test_compile_resolved.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr
    env=s.V.environment();env['VELA_MERIT_FIXTURE_MANIFEST']=str(LOCAL/'fixtures.tsv');r=subprocess.run([str(LOCAL/'test.exe')],env=env,capture_output=True,text=True);(LOCAL/'test_result.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stdout+r.stderr
    seal();print(r.stdout,flush=True)

def run():
    s.a.verify(OUT/'build_evidence.json');s.a.verify(w.OUT/'settled_source/dc_evidence.json');jobs=[];files=[OUT/'build_evidence.json',w.OUT/'settled_source/dc_evidence.json']
    controls=[('failed_e1089',s.LOCAL/'restart_coordinates/dc/native_reload'),('qualified_control',s.prior.LOCAL/'candidate/dc/phumob/m65_n23_vd_1p000000_endpoint/vela/vg_010/attempt_0')]
    for label,base in controls:
        dest=LOCAL/'repair'/label;cfg=s.a.read(base/'config.json');cfg.update(state_file=str(base/'state.csv'),output_state_file=str(dest/'state.csv'));cfg['solver']['local_update_diagnostics']=dict(enabled=True,nodes=[967,983,1089],csv_file=str(dest/'updates.csv'),first_iterations=200,every_iterations=1)
        s.a.write(dest/'config.json',cfg);s.V.post_config(cfg,dest);files += [base/'state.csv']+list(dest.glob('*.json'));jobs.append(dict(key=label,device='n23',dest=str(dest),kind='repair'))
    for c in s.a.read(w.OUT/'settled_source/contract.json')['cases']:
        for j in c['jobs']:
            dest=LOCAL/'source'/c['key']/j['label'];cfg=s.a.read(Path(j['dest'])/'config.json');cfg['output_state_file']=str(dest/'state.csv');s.a.write(dest/'config.json',cfg);s.V.post_config(cfg,dest);files += [Path(cfg['state_file'])]+list(dest.glob('*.json'));jobs.append(dict(key=c['key'],device=c['device'],dest=str(dest),kind='source',source_file=c['source_file'],alpha=j['alpha'],label=j['label']))
    s.a.write(OUT/'contract.json',dict(jobs=jobs,change='Only replace unweighted L2 norm decrease comparison with exact-sign squared-norm difference. Existing zero-to-zero exception, carrier validity, damping grid, refinements, caps and row gate unchanged.',frozen_row_gate=1e-6,kcl_over_Id=1e-8,production_modified=False,finite_volume_replacement=False))
    s.d.matrix.freeze(OUT/'freeze.json',files+[OUT/'contract.json'])
    def one(j):
        dest=Path(j['dest']);env=s.env(j,j['alpha']) if j['kind']=='source' else s.V.environment();env.update(VELA_STABLE_MERIT='1',VELA_STABLE_MERIT_TRACE=str(dest/'merit_trace.csv'))
        status=s.V.execute(dest/'config.json',RUNNER,env)
        for name in ('all_row','acceptance_edges'):s.V.execute(dest/(name+'.json'),RUNNER,env)
        result=dict(**j,**s.prior.q.run.w.old.prior.old.qualify(j,dest));print(result,flush=True);return result
    with ThreadPoolExecutor(max_workers=2) as pool:results=list(pool.map(one,jobs))
    columns=list(dict.fromkeys(k for r in results for k in r));s.a.write_csv(OUT/'dc.csv',[{k:r.get(k,'') for k in columns} for r in results]);s.d.matrix.freeze(OUT/'dc_evidence.json',[OUT/'freeze.json',OUT/'dc.csv']+[p for j in jobs for p in Path(j['dest']).rglob('*') if p.is_file()])

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('build','finish','run'));globals()[p.parse_args().action]()
