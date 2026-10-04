"""Output-only current-assembler/current-Newton trial operands at electron 1089."""
import argparse
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import calibrate_simplemos_phumob_local_source_20260909 as s
import audit_simplemos_linesearch_precision_20260908 as old

LOCAL=s.REPO/'build-release/phumob_restart_followup_20260909/linesearch'
OUT=s.REPO/'reference_tcad/simplemos_sentaurus2022/phumob_restart_followup_20260909/linesearch'
RUNNER=LOCAL/'runner.exe'

def build():
    s.a.verify(s.OUT/'restart_coordinates/build_evidence.json');LOCAL.mkdir(parents=True,exist_ok=False)
    (LOCAL/'audit_shared.h').write_text('#pragma once\n#include <fstream>\n#include <iomanip>\n#include <string>\n#include <cstdlib>\nnamespace vela_lsa { extern std::string prefix; }\n')
    source=(s.REPO/'src/equation/CoupledDDAssembler.cpp').read_text(encoding='utf-8')
    start=source.index('VectorXd CoupledDDAssembler::residualImpl(');end=source.index('std::vector<CoupledDDCarrierTermDiagnostic>',start);frag=source[start:end];assert frag.count('    return r;')==1
    source=source[:start]+frag.replace('    return r;',old.POISSON_DUMP+'\n    return r;')+source[end:]
    start=source.index('std::vector<CoupledDDEdgeFluxDiagnostic>\nCoupledDDAssembler::sgEdgeFluxDiagnostics(');end=source.index('std::vector<CoupledDDElectronTransportFactorDiagnostic>',start);frag=source[start:end]
    node=old.kernel.NODE_DUMP.replace('if(const char* auditPrefix=std::getenv("VELA_MINORITY_PRECISION_PREFIX")) {','if(!vela_lsa::prefix.empty()) { const char* auditPrefix=vela_lsa::prefix.c_str();')
    marker='    std::vector<CoupledDDEdgeFluxDiagnostic> edges;';assert frag.count(marker)==1;frag=frag.replace(marker,node+'\n'+marker)
    assert frag.count('        edges.push_back(record);')==1;frag=frag.replace('        edges.push_back(record);',old.kernel.EDGE_DUMP+'\n        edges.push_back(record);')
    (LOCAL/'CoupledDDAssembler.cpp').write_text('#include "audit_shared.h"\n'+source[:start]+frag+source[end:],encoding='utf-8')
    source=(s.LOCAL/'restart_coordinates/NewtonSolver.cpp').read_text(encoding='utf-8')
    start=source.index('        const auto runLineSearch = [&](const VectorXd& trialStep) {');end=source.index('        auto ls = runLineSearch(step);',start);frag=source[start:end]
    assert frag.count('            return lineSearch.search(')==1;frag=frag.replace('            return lineSearch.search(',old.SEARCH_PREFIX+'\n            return lineSearch.search(')
    assert frag.count('                    return assembler.residual(candidate, bcs);')==1;frag=frag.replace('                    return assembler.residual(candidate, bcs);',old.SEARCH_RESIDUAL)
    (LOCAL/'NewtonSolver.cpp').write_text('#include "audit_shared.h"\nnamespace vela_lsa { std::string prefix; }\n'+source[:start]+frag+source[end:],encoding='utf-8')
    commands=[]
    for name,path in (('CoupledDDAssembler',s.LOCAL/'scaled_source/compile_command.json'),('NewtonSolver',s.LOCAL/'restart_coordinates/compile_command.json')):
        cmd=s.a.read(path);cmd[cmd.index('-c')+1]=str(LOCAL/(name+'.cpp'));cmd[cmd.index('-o')+1]=str(LOCAL/(name+'.o'));commands.append(cmd)
    s.a.write(LOCAL/'compile_commands.json',commands)
    def compile_one(cmd):
        r=subprocess.run(cmd,cwd=s.REPO/'build-release',env=s.V.environment(),capture_output=True,text=True);(LOCAL/(Path(cmd[cmd.index('-c')+1]).stem+'.compile.log')).write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr[-4000:]
    with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(compile_one,commands))
    cmd=s.a.read(s.LOCAL/'scaled_source/link_command.json');cmd[1]=str(LOCAL/'CoupledDDAssembler.o');cmd.insert(1,str(LOCAL/'NewtonSolver.o'));cmd[cmd.index('-o')+1]=str(RUNNER)
    s.a.write(LOCAL/'link_command.json',cmd);r=subprocess.run(cmd,env=s.V.environment(),capture_output=True,text=True);(LOCAL/'link.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr
    s.d.matrix.freeze(OUT/'build_evidence.json',[Path(__file__).resolve(),Path(old.__file__),Path(old.kernel.__file__),s.REPO/'src/equation/CoupledDDAssembler.cpp',s.OUT/'restart_coordinates/build_evidence.json',s.REPO/'build-release/libvela_core.a']+[p for p in LOCAL.iterdir() if p.is_file()])

def run():
    s.a.verify(OUT/'build_evidence.json');jobs=[];files=[OUT/'build_evidence.json']
    bases=[('failed_e1089',s.LOCAL/'restart_coordinates/dc/native_reload'),('qualified_control',s.prior.LOCAL/'candidate/dc/phumob/m65_n23_vd_1p000000_endpoint/vela/vg_010/attempt_0')]
    for label,base in bases:
        for mode in ('plain','trace'):
            dest=LOCAL/label/mode;cfg=s.a.read(base/'config.json');cfg.update(state_file=str(base/'state.csv'),output_state_file=str(dest/'state.csv'))
            cfg['solver']['max_iter']=1;cfg['solver']['carrier_row_convergence']['min_newton_max_iter']=1
            cfg['solver']['local_update_diagnostics']=dict(enabled=True,nodes=[967,983,1089],csv_file=str(dest/'updates.csv'),first_iterations=1,every_iterations=1)
            s.a.write(dest/'config.json',cfg);s.V.post_config(cfg,dest);files += [base/'state.csv',base/'config.status.json']+list(dest.glob('*.json'));jobs.append(dict(label=label,mode=mode,base=str(base),dest=str(dest)))
    s.a.write(OUT/'contract.json',dict(jobs=jobs,scope='One actual Newton iteration from final e1089 failed state plus same-bias qualified control; trace on/off identity. Preserve four refinements, step caps, production line search. Diagnostic iteration budget is not a new qualification attempt.',node=1089,original_row_gate=1e-6,production_changed=False))
    s.d.matrix.freeze(OUT/'freeze.json',files+[OUT/'contract.json']);rows=[]
    for j in jobs:
        dest=Path(j['dest']);env=s.V.environment()
        if j['mode']=='trace':env.update(VELA_LS_AUDIT_DIR=str(dest),VELA_LS_AUDIT_ITERATION='1')
        status=s.V.execute(dest/'config.json',RUNNER,env)
        for name in ('all_row','acceptance_edges'):s.V.execute(dest/(name+'.json'),RUNNER,s.V.environment())
        rows.append(dict(label=j['label'],mode=j['mode'],**s.prior.q.run.w.old.prior.old.qualify(dict(device='n23'),dest)));print(rows[-1],flush=True)
    identities=[]
    for label,_ in bases:
        plain=LOCAL/label/'plain';trace=LOCAL/label/'trace';p=s.a.read(plain/'config.status.json');t=s.a.read(trace/'config.status.json')
        identical=s.a.sha(plain/'state.csv')==s.a.sha(trace/'state.csv') and all(p[k]==t[k] for k in ('iterations','exit_code','failure_reason','converged','final_residual'))
        dumps=list(trace.glob('*_direction.csv'));repeat=all(all(r['R']==r['repeated_R'] for r in s.a.rows(f)) for f in dumps)
        identities.append(dict(label=label,state_and_status_identity=identical,residual_repeat_identity=repeat,snapshots=len(dumps),qualified=identical and repeat));assert identical and repeat
    s.a.write_csv(OUT/'dc.csv',rows);s.a.write_csv(OUT/'identity.csv',identities);s.d.matrix.freeze(OUT/'replay_evidence.json',[OUT/'freeze.json',OUT/'dc.csv',OUT/'identity.csv']+[p for label,_ in bases for p in (LOCAL/label).rglob('*') if p.is_file()])

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('build','run'));globals()[p.parse_args().action]()
