"""Build a guarded merit-comparison-only Newton overlay and property tests."""
from pathlib import Path
import subprocess
import audit_simplemos_linesearch_precision_20260908 as audit
import run_simplemos_linesearch_precision_20260908 as precision

t=audit.t;c=t.c;a=t.a;d=t.d;v=t.v;p=t.p
LOCAL=p.REPO/'build-release/simplemos_stable_merit_20260908'
OUT=p.REPO/'reference_tcad/simplemos_sentaurus2022/stable_merit_20260908'
RUNNER=LOCAL/'runner.exe'
HEADER=p.REPO/'scripts/diagnostics/simplemos_stable_merit.hpp'
TEST=p.REPO/'tests/diagnostics/test_simplemos_stable_merit.cpp'
PREFIX=r'''
        const char* stableText=std::getenv("VELA_STABLE_MERIT");
        const bool stableMerit=stableText && std::string(stableText)=="1";
        if(stableText && !stableMerit)throw std::runtime_error("Invalid stable merit selector");
        if(stableMerit && (cfg_.lineSearchMode!="merit" || cfg_.globalContinuityClosure.mode!="off" ||
            (cfg_.residualNorm!="l2" && (residualScales.psi!=1. || residualScales.phin!=1. || residualScales.phip!=1. ||
              residualWeights.psi!=1. || residualWeights.phin!=1. || residualWeights.phip!=1.))))
            throw std::runtime_error("Stable merit audit requires unweighted L2 and no global merit term");
        const auto stableDecrease=[&](const VectorXd& candidateResidual,Real alpha) {
            const auto check=simplemos_stable_merit::compare(r,candidateResidual);
            if(const char* name=std::getenv("VELA_STABLE_MERIT_TRACE")) {
                const bool header=!std::filesystem::exists(name);std::ofstream f(name,std::ios::app);
                if(header)f<<"iteration,alpha,accepted,sign,exact_fallback,delta,roundoff_bound\n";
                f<<std::setprecision(21)<<iter<<','<<alpha<<','<<check.accepted<<','<<check.sign<<','<<check.exactFallback<<','<<check.delta<<','<<check.roundoffBound<<'\n';
            }
            return check.accepted;
        };
'''


def main():
    a.verify(audit.OUT/'final_evidence.json');a.verify(t.OUT/'build_evidence.json')
    LOCAL.mkdir(parents=True,exist_ok=True);OUT.mkdir(parents=True,exist_ok=True);assert not (OUT/'build_evidence.json').exists()
    if RUNNER.exists():
        finish_tests(a.read(LOCAL/'compile_command.json'),a.read(LOCAL/'link_command.json'));return
    cmds=a.read(c.old.LOCAL/'compile_commands.json');cmd=list(next(x for x in cmds if Path(x[x.index('-c')+1]).name=='NewtonSolver.cpp'))
    oldobj=cmd[cmd.index('-o')+1];source=Path(cmd[cmd.index('-c')+1]).read_text(encoding='utf-8')
    marker='        const auto runLineSearch = [&](const VectorXd& trialStep) {';assert source.count(marker)==1;source=source.replace(marker,PREFIX+'\n'+marker)
    start=source.index(marker);end=source.index('        auto ls = runLineSearch(step);',start);frag=source[start:end]
    target='                    : BacktrackingLineSearch::DecreaseAcceptFunction{});';assert frag.count(target)==1
    frag=frag.replace(target,'                    : (stableMerit ? BacktrackingLineSearch::DecreaseAcceptFunction(stableDecrease) : BacktrackingLineSearch::DecreaseAcceptFunction{}));')
    source=source[:start]+frag+source[end:];(LOCAL/'NewtonSolver.cpp').write_text('#include "simplemos_stable_merit.hpp"\n#include <cstdlib>\n'+source,encoding='utf-8')
    cmd.insert(1,'-I'+str(HEADER.parent));cmd[cmd.index('-c')+1]=str(LOCAL/'NewtonSolver.cpp');cmd[cmd.index('-o')+1]=str(LOCAL/'NewtonSolver.o')
    a.write(LOCAL/'compile_command.json',cmd)
    result=subprocess.run(cmd,cwd=p.REPO/'build-release',env=v.environment(),capture_output=True,encoding='utf-8',errors='replace');(LOCAL/'build.log').write_text(result.stdout+result.stderr,encoding='utf-8');assert result.returncode==0,result.stderr[-4000:]
    link=[str(LOCAL/'NewtonSolver.o') if x==oldobj else x for x in a.read(t.LOCAL/'link_command.json')];link[link.index('-o')+1]=str(RUNNER)
    a.write(LOCAL/'link_command.json',link);result=subprocess.run(link,env=v.environment(),capture_output=True,encoding='utf-8',errors='replace');(LOCAL/'link.log').write_text(result.stdout+result.stderr);assert result.returncode==0,result.stderr[-3000:]
    finish_tests(cmd,link)


def finish_tests(cmd,link):
    fixtures=[]
    for job in a.read(audit.OUT/'contract.json')['jobs']:
        root=Path(job['dest'])
        for i in range(13):
            label='trial_'+str(i);r=next(x for x in a.rows(precision.core.OUT/'merit.csv') if x['key']==job['case']['key'] and x['mode']=='production' and x['label']==label)
            fixtures.append(str(root/'base_direction.csv')+'\t'+str(root/(label+'_direction.csv'))+'\t'+('1' if r['would_decrease']=='True' else '0'))
    (LOCAL/'fixtures.tsv').write_text('\n'.join(fixtures)+'\n')
    flags=[x for x in cmd[1:] if x.startswith(('-I','-D','-O','-std='))];testcmd=[cmd[0]]+flags+['-ID:/msys64/ucrt64/include/eigen3',str(TEST),'D:/msys64/ucrt64/lib/libCatch2Main.a','D:/msys64/ucrt64/lib/libCatch2.a','-o',str(LOCAL/'test.exe')]
    a.write(LOCAL/'test_compile_resolved.json',testcmd);result=subprocess.run(testcmd,env=v.environment(),capture_output=True,encoding='utf-8',errors='replace');(LOCAL/'test_build_resolved.log').write_text(result.stdout+result.stderr,encoding='utf-8');assert result.returncode==0,result.stderr[-3000:]
    e=v.environment();e['VELA_MERIT_FIXTURE_MANIFEST']=str(LOCAL/'fixtures.tsv');result=subprocess.run([str(LOCAL/'test.exe')],env=e,capture_output=True,encoding='utf-8',errors='replace');(LOCAL/'test_result.log').write_text(result.stdout+result.stderr);assert result.returncode==0,result.stdout+result.stderr
    d.matrix.freeze(OUT/'build_evidence.json',[Path(__file__).resolve(),HEADER,TEST,audit.OUT/'final_evidence.json',t.OUT/'build_evidence.json']+[x for x in LOCAL.rglob('*') if x.is_file()]+[Path(x) for x in link if x.endswith('.o')]);print(result.stdout,flush=True)


if __name__=='__main__':main()
