"""Isolate higher-precision Poisson evaluation without changing state storage."""
import argparse
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import audit_simplemos_lowvd_source_precision_20260910 as v
s=v.s;m=v.m
LOCAL=v.LOCAL/'poisson_candidate';OUT=v.OUT/'poisson_candidate';RUNNER=LOCAL/'runner.exe'
HEADER=s.REPO/'scripts/diagnostics/simplemos_poisson_precision.hpp'
HOOK=r'''
    if (const char* precisionText=std::getenv("VELA_POISSON_PRECISION")) {
        const std::string precisionMode(precisionText);
        if(precisionMode!="kernel" && precisionMode!="packed")throw std::runtime_error("Invalid Poisson precision mode");
        if(usesFermiDirac_ || substitution!=nullptr || !bcs.thermionic.empty())throw std::runtime_error("Unsupported Poisson precision diagnostic branch");
        namespace pp=simplemos_poisson_precision;
        std::vector<pp::Node> pn;pn.reserve(N);
        for(int i=0;i<N;++i)pn.push_back({psi(i),n(i),p(i),x(i),x(N+i),x(2*N+i),potentialScale,
            electronQuasiFermiReferenceAt(i),holeQuasiFermiReferenceAt(i),ni_[i],Vt_,
            doping_.netDoping(i),poissonElectronVol_[i],poissonHoleVol_[i],poissonDopantVol_[i],
            constants::q,chargeAreaFactor,scaling_.enabled?scaling_.permittivityReference_F_per_m*scaling_.V0:1.,fixedInterfaceChargeRhs_(i)});
        std::vector<pp::Edge> pe;pe.reserve(mesh_.numEdges());
        for(Index e=0;e<mesh_.numEdges();++e) {
            const auto& edge=mesh_.getEdge(e);if(edge.length<1e-30)continue;
            pe.push_back({static_cast<int>(edge.n0),static_cast<int>(edge.n1),edgeAssemblyKernels_[e].poissonCoupling});
        }
        const auto pr=pp::evaluate(pn,pe,precisionMode=="packed");
        for(int i=0;i<N;++i)r(i)=static_cast<Real>(pr[i]);
    }
'''

def build():
    s.a.verify(v.OUT/'replay_evidence.json');LOCAL.mkdir(parents=True,exist_ok=False)
    source=(s.LOCAL/'scaled_source/CoupledDDAssembler.cpp').read_text(encoding='utf-8')
    start=source.index('VectorXd CoupledDDAssembler::residualImpl(');end=source.index('std::vector<CoupledDDCarrierTermDiagnostic>',start);frag=source[start:end]
    marker='    for (const auto& [node, value] : bcs.psi)';assert frag.count(marker)==1
    source='#include "simplemos_poisson_precision.hpp"\n'+source[:start]+frag.replace(marker,HOOK+'\n'+marker)+source[end:]
    (LOCAL/'CoupledDDAssembler.cpp').write_text(source,encoding='utf-8')
    cmd=s.a.read(s.LOCAL/'scaled_source/compile_command.json');cmd.insert(1,'-I'+str(HEADER.parent));cmd[cmd.index('-c')+1]=str(LOCAL/'CoupledDDAssembler.cpp');cmd[cmd.index('-o')+1]=str(LOCAL/'assembler.o');s.a.write(LOCAL/'compile_command.json',cmd)
    r=subprocess.run(cmd,cwd=s.REPO/'build-release',env=s.V.environment(),capture_output=True,text=True);(LOCAL/'compile.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr[-4000:]
    link=s.a.read(m.LOCAL/'link_command.json');old=str(s.LOCAL/'scaled_source/assembler.o');assert link.count(old)==1;link[link.index(old)]=str(LOCAL/'assembler.o');link[link.index('-o')+1]=str(RUNNER);s.a.write(LOCAL/'link_command.json',link)
    r=subprocess.run(link,env=s.V.environment(),capture_output=True,text=True);(LOCAL/'link.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr[-4000:]
    s.d.matrix.freeze(OUT/'build_evidence.json',[Path(__file__).resolve(),HEADER,m.OUT/'build_evidence.json',v.OUT/'replay_evidence.json']+[p for p in LOCAL.iterdir() if p.is_file()])

def run():
    s.a.verify(OUT/'build_evidence.json');s.a.verify(v.OUT/'precision/evidence.json')
    jobs=[];files=[OUT/'build_evidence.json',v.OUT/'precision/evidence.json']
    for case in s.a.read(m.OUT/'contract.json')['jobs']:
        if case['kind']!='source' or case['label'] not in ('zero','plus_full','minus_full'):continue
        for mode in ('kernel','packed'):
            base=Path(case['dest']);dest=LOCAL/mode/case['key']/case['label'];cfg=s.a.read(base/'config.json');cfg['output_state_file']=str(dest/'state.csv')
            assert cfg['solver'].get('carrier_statistics',{}).get('model','boltzmann')=='boltzmann'
            s.a.write(dest/'config.json',cfg);s.V.post_config(cfg,dest)
            jobs.append(dict(case,base=str(base),dest=str(dest),mode=mode));files += list(dest.glob('*.json'))+[Path(cfg['state_file'])]
    s.a.write(OUT/'contract.json',dict(jobs=jobs,change='Replace only Poisson residual evaluation with binary128 arithmetic. Kernel uses the existing double psi/n/p; packed recomputes physical psi and Boltzmann density from the actual packed double state. Original Jacobian, state storage, stable merit, source, 200 iterations, 1e-6 carrier gate and 1e-8 KCL gate retained.',production_changed=False,qualification='A diagnostic candidate. Do not combine qualified points with prior implementations to claim 16/16.'))
    s.d.matrix.freeze(OUT/'freeze.json',files+[OUT/'contract.json'])
    def one(j):
        dest=Path(j['dest']);env=s.env(j,j['alpha']);env.update(VELA_STABLE_MERIT='1',VELA_POISSON_PRECISION=j['mode'])
        s.V.execute(dest/'config.json',RUNNER,env)
        for name in ('all_row','acceptance_edges'):s.V.execute(dest/(name+'.json'),RUNNER,env)
        result=dict(**j,**s.prior.q.run.w.old.prior.old.qualify(j,dest));print(result,flush=True);return result
    with ThreadPoolExecutor(max_workers=2) as pool:results=list(pool.map(one,jobs))
    s.a.write_csv(OUT/'dc.csv',results)
    s.d.matrix.freeze(OUT/'dc_evidence.json',[OUT/'freeze.json',OUT/'dc.csv']+[p for j in jobs for p in Path(j['dest']).rglob('*') if p.is_file()])

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('build','run'));globals()[p.parse_args().action]()
