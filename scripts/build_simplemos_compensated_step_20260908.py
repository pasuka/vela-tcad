"""Second isolated numerical candidate: accumulate lost psi updates only."""
from pathlib import Path
import subprocess
from concurrent.futures import ThreadPoolExecutor
import build_simplemos_stable_merit_20260908 as b

t=b.t;c=b.c;a=b.a;d=b.d;v=b.v;p=b.p
LOCAL=p.REPO/'build-release/simplemos_compensated_step_20260908'
OUT=p.REPO/'reference_tcad/simplemos_sentaurus2022/compensated_step_20260908'
RUNNER=LOCAL/'runner.exe';HEADER=p.REPO/'scripts/diagnostics/simplemos_compensated_step.hpp';TEST=p.REPO/'tests/diagnostics/test_simplemos_compensated_step.cpp'


def main():
    a.verify(b.OUT/'build_evidence.json');LOCAL.mkdir(parents=True,exist_ok=False);OUT.mkdir(parents=True,exist_ok=False)
    ns=(b.LOCAL/'NewtonSolver.cpp').read_text(encoding='utf-8');marker='    LinearSolver linearSolver;\n    LineSearchConfig lscfg;';assert ns.count(marker)==1
    init='''    const char* carryText=std::getenv("VELA_COMPENSATED_PSI_STEP");
    if(carryText && std::string(carryText)!="1")throw std::runtime_error("Invalid compensated psi selector");
    simplemos_compensated_step::active=carryText!=nullptr;
    simplemos_compensated_step::carry.assign(N,0.);
'''
    ns='#include "simplemos_compensated_step.hpp"\nnamespace simplemos_compensated_step { bool active=false; std::vector<double> carry; }\n'+ns.replace(marker,init+marker)
    (LOCAL/'NewtonSolver.cpp').write_text(ns,encoding='utf-8')
    ls=(p.REPO/'src/numerics/LineSearch.cpp').read_text(encoding='utf-8');marker='        VectorXd candidate = x + alpha * step;';assert ls.count(marker)==1
    ls=ls.replace(marker,marker+'''
        std::vector<double> trialCarry;
        if(simplemos_compensated_step::active) {
            trialCarry.resize(simplemos_compensated_step::carry.size());
            if(x.size()!=3*static_cast<int>(trialCarry.size()))throw std::runtime_error("Compensated psi state size mismatch");
            for(std::size_t i=0;i<trialCarry.size();++i) {
                const auto p=simplemos_compensated_step::propose(x(i),step(i),alpha,simplemos_compensated_step::carry[i]);
                candidate(i)=p.value;trialCarry[i]=p.remainder;
            }
        }
''')
    marker='        if (accepted) {\n            incrementPerformanceCounter("newton.line_search_accepted");';assert ls.count(marker)==1
    ls=ls.replace(marker,'''        if (accepted) {
            if(simplemos_compensated_step::active)simplemos_compensated_step::carry=std::move(trialCarry);
            incrementPerformanceCounter("newton.line_search_accepted");''')
    (LOCAL/'LineSearch.cpp').write_text('#include "simplemos_compensated_step.hpp"\n'+ls,encoding='utf-8')
    oldcmds=a.read(c.old.LOCAL/'compile_commands.json');cmds=[list(a.read(b.LOCAL/'compile_command.json')),list(next(x for x in oldcmds if Path(x[x.index('-c')+1]).name=='LineSearch.cpp'))];replace={}
    for cmd in cmds:
        stem=Path(cmd[cmd.index('-c')+1]).stem;replace[cmd[cmd.index('-o')+1]]=str(LOCAL/(stem+'.o'));cmd.insert(1,'-I'+str(HEADER.parent));cmd[cmd.index('-c')+1]=str(LOCAL/(stem+'.cpp'));cmd[cmd.index('-o')+1]=str(LOCAL/(stem+'.o'))
    a.write(LOCAL/'compile_commands.json',cmds)
    def compile(cmd):
        r=subprocess.run(cmd,cwd=p.REPO/'build-release',env=v.environment(),capture_output=True,encoding='utf-8',errors='replace');(LOCAL/(Path(cmd[cmd.index('-c')+1]).stem+'.build.log')).write_text(r.stdout+r.stderr,encoding='utf-8');assert r.returncode==0,r.stderr[-3000:]
    with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(compile,cmds))
    link=[replace.get(x,x) for x in a.read(b.LOCAL/'link_command.json')];link[link.index('-o')+1]=str(RUNNER);a.write(LOCAL/'link_command.json',link)
    r=subprocess.run(link,env=v.environment(),capture_output=True,encoding='utf-8',errors='replace');(LOCAL/'link.log').write_text(r.stdout+r.stderr,encoding='utf-8');assert r.returncode==0,r.stderr[-3000:]
    flags=[x for x in cmds[0] if x.startswith(('-I','-D','-O','-std='))];cmd=[cmds[0][0]]+flags+['-ID:/msys64/ucrt64/include/eigen3',str(TEST),'D:/msys64/ucrt64/lib/libCatch2Main.a','D:/msys64/ucrt64/lib/libCatch2.a','-o',str(LOCAL/'test.exe')]
    a.write(LOCAL/'test_compile.json',cmd);r=subprocess.run(cmd,env=v.environment(),capture_output=True,encoding='utf-8',errors='replace');(LOCAL/'test_build.log').write_text(r.stdout+r.stderr,encoding='utf-8');assert r.returncode==0,r.stderr[-3000:]
    r=subprocess.run([str(LOCAL/'test.exe')],env=v.environment(),capture_output=True,encoding='utf-8',errors='replace');(LOCAL/'test_result.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stdout+r.stderr
    d.matrix.freeze(OUT/'build_evidence.json',[Path(__file__).resolve(),HEADER,TEST,b.OUT/'build_evidence.json']+[x for x in LOCAL.rglob('*') if x.is_file()]+[Path(x) for x in link if x.endswith('.o')]);print(r.stdout,flush=True)


if __name__=='__main__':main()
