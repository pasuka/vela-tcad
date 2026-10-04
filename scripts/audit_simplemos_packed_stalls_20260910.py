"""Output-only trial replay after exact packed checkpoint preservation."""
import argparse
import shlex
import subprocess
from pathlib import Path
import validate_simplemos_packed_restart_v2_20260910 as v
a,d=v.a,v.d
LOCAL=v.LOCAL/'stalls';OUT=v.OUT/'stalls';RUNNER=LOCAL/'runner.exe'
HOOK=r'''
                    const VectorXd result=assembler.residual(candidate,bcs);
                    if(const char* directory=std::getenv("VELA_PACK_STALL")) {
                        static int sequence=0;
                        const std::string prefix=std::string(directory)+"/trial_"+std::to_string(sequence++);
                        std::ofstream out(prefix+".csv");out<<std::setprecision(17);
                        out<<"row,x,raw_step,capped_step,candidate,R,Rtrial,JdxR,scale,weight\n";
                        const Real scales[3]={residualScales.psi,residualScales.phin,residualScales.phip};
                        const Real weights[3]={residualWeights.psi,residualWeights.phin,residualWeights.phip};
                        for(int i=0;i<3*N;++i)out<<i<<','<<x(i)<<','<<rawStep(i)<<','<<trialStep(i)<<','<<candidate(i)<<','<<r(i)<<','<<result(i)<<','<<rawLinearResidual(i)<<','<<scales[i/N]<<','<<weights[i/N]<<'\n';
                    }
                    return result;
'''

def build():
    a.verify(v.OUT/'cohort_v2/validation_freeze.json');LOCAL.mkdir(parents=True,exist_ok=False)
    src=v.REPO/'src/solver/NewtonSolver.cpp';text=src.read_text(encoding='utf-8')
    start=text.index('        const auto runLineSearch = [&](const VectorXd& trialStep) {')
    end=text.index('        auto ls = runLineSearch(step);',start);frag=text[start:end]
    marker='                    return assembler.residual(candidate, bcs);';assert frag.count(marker)==1
    target=LOCAL/'NewtonSolver.cpp';target.write_text('#include <cstdlib>\n'+text[:start]+frag.replace(marker,HOOK)+text[end:],encoding='utf-8')
    entry=next(e for e in a.read(v.REPO/'build-release/compile_commands.json') if Path(e['file'])==src)
    cmd=shlex.split(entry['command'].replace('\\','/'));cmd=[('-DVELA_VERSION="0.1.0"' if x.startswith('-DVELA_VERSION=') else x) for x in cmd]
    cmd[cmd.index('-c')+1]=str(target);cmd[cmd.index('-o')+1]=str(LOCAL/'newton.o');a.write(LOCAL/'compile.json',cmd)
    p=subprocess.run(cmd,cwd=entry['directory'],env=v.run.V.environment(),capture_output=True,text=True);(LOCAL/'compile.log').write_text(p.stdout+p.stderr);assert p.returncode==0,p.stderr[-3000:]
    link=[cmd[0],str(LOCAL/'newton.o'),str(v.REPO/'build-release/CMakeFiles/vela_example_runner.dir/src/tools/vela_example_runner.cpp.obj'),str(v.REPO/'build-release/libvela_core.a')]
    link += ['D:/msys64/ucrt64/lib/lib'+n for n in ('spdlog.dll.a','fmt.a','umfpack.dll.a','spqr.dll.a','cholmod.dll.a')]+['-o',str(RUNNER)]
    a.write(LOCAL/'link.json',link);p=subprocess.run(link,env=v.run.V.environment(),capture_output=True,text=True);(LOCAL/'link.log').write_text(p.stdout+p.stderr);assert p.returncode==0,p.stderr
    d.matrix.freeze(OUT/'build_evidence.json',[Path(__file__).resolve(),src,v.OUT/'cohort_v2/validation_freeze.json',v.REPO/'build-release/libvela_core.a',v.REPO/'build-release/CMakeFiles/vela_example_runner.dir/src/tools/vela_example_runner.cpp.obj']+[p for p in LOCAL.iterdir() if p.is_file()])

def replay():
    a.verify(OUT/'build_evidence.json');a.verify(v.OUT/'cohort_v2/comparison_evidence.json')
    files=[OUT/'build_evidence.json',v.OUT/'cohort_v2/comparison_evidence.json'];rows=[]
    for j in a.rows(v.OUT/'cohort_v2/selected_states.csv'):
        if j['qualified']=='True':continue
        base=Path(j['dest']);dest=LOCAL/j['device']/str(j['vd'])/str(j['index'])
        cfg=a.read(base/'config.json');cfg['output_state_file']=str(dest/'state.csv')
        cfg['solver']['local_update_diagnostics']=dict(enabled=True,nodes=[int(a.read(base/'config.status.json')['carrier_row_convergence']['max_ratio_node'])],csv_file=str(dest/'updates.csv'),first_iterations=1,every_iterations=1)
        a.write(dest/'config.json',cfg);env=v.run.V.environment();env['VELA_PACK_STALL']=str(dest)
        status=v.run.V.execute(dest/'config.json',RUNNER,env);old=a.read(base/'config.status.json')
        same=a.sha(base/'state.csv')==a.sha(dest/'state.csv') and all(old[k]==status[k] for k in ('iterations','final_residual','failure_reason','exit_code'))
        rows.append(dict(device=j['device'],vd=j['vd'],index=j['index'],identity=same,trials=len(list(dest.glob('trial_*.csv'))),dest=str(dest)))
        assert same,rows[-1];print(rows[-1],flush=True)
        files += [base/'config.json',base/'state.csv',base/'config.status.json']+[p for p in dest.iterdir() if p.is_file()]
    a.write_csv(OUT/'identity.csv',rows);d.matrix.freeze(OUT/'replay_evidence.json',files+[OUT/'identity.csv'])

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('build','replay'));globals()[p.parse_args().action]()
