"""Output-only identity replays of the final rejected Newton directions."""
import argparse
import copy
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import restart_simplemos_srh_finite_20260908 as prior
import audit_simplemos_minority_residual_20260906 as kernel

t=prior.t; c=t.c; a=t.a; d=t.d; v=t.v; p=t.p
LOCAL=p.REPO/'build-release/simplemos_linesearch_srh_audit_20260908'
OUT=p.REPO/'reference_tcad/simplemos_sentaurus2022/linesearch_srh_audit_20260908'
RUNNER=LOCAL/'runner.exe'

POISSON_DUMP=r'''
    if (!vela_lsa::prefix.empty()) {
        std::ofstream ed(vela_lsa::prefix+"_poisson_edges.csv");
        ed<<std::setprecision(17)<<"node0,node1,G,rounded_flux\n";
        for(Index e=0;e<mesh_.numEdges();++e) {
            const auto& edge=mesh_.getEdge(e);if(edge.length<1e-30)continue;
            const Real G=edgeAssemblyKernels_[e].poissonCoupling;
            ed<<edge.n0<<','<<edge.n1<<','<<G<<','<<G*(psi(edge.n0)-psi(edge.n1))<<'\n';
        }
        std::ofstream nd(vela_lsa::prefix+"_poisson_nodes.csv");
        nd<<std::setprecision(17)<<"node,xpsi,xn,xp,psi,n,p,net_doping,voln,volp,vold,q,area_factor,poisson_scale,potential_scale,interface_rhs,bcpsi,bcn,bcp,Rpsi,Rn,Rp,rounded_hole_charge,rounded_dopant_charge,rounded_electron_charge\n";
        for(int i=0;i<N;++i)nd<<i<<','<<x(i)<<','<<x(N+i)<<','<<x(2*N+i)<<','<<psi(i)<<','<<n(i)<<','<<p(i)<<','<<doping_.netDoping(i)<<','
          <<poissonElectronVol_[i]<<','<<poissonHoleVol_[i]<<','<<poissonDopantVol_[i]<<','<<constants::q<<','<<chargeAreaFactor<<','
          <<(scaling_.enabled?scaling_.permittivityReference_F_per_m*scaling_.V0:1.)<<','<<potentialScale<<','<<fixedInterfaceChargeRhs_(i)<<','
          <<bcs.psi.count(i)<<','<<bcs.phin.count(i)<<','<<bcs.phip.count(i)<<','<<r(i)<<','<<r(N+i)<<','<<r(2*N+i)<<','
          <<constants::q*p(i)*poissonHoleVol_[i]*chargeAreaFactor<<','<<constants::q*doping_.netDoping(i)*poissonDopantVol_[i]*chargeAreaFactor<<','
          <<constants::q*n(i)*poissonElectronVol_[i]*chargeAreaFactor<<'\n';
    }
'''

SEARCH_PREFIX=r'''
            const char* auditDir=std::getenv("VELA_LS_AUDIT_DIR");
            const char* auditIter=std::getenv("VELA_LS_AUDIT_ITERATION");
            const bool auditEnabled=auditDir && auditIter && iter==std::stoi(auditIter);
            int auditAttempt=0;
            auto auditDump=[&](const VectorXd& state,const VectorXd& observed,const std::string& label,Real alpha) {
                vela_lsa::prefix=std::string(auditDir)+"/"+label;
                const VectorXd repeated=assembler.residual(state,bcs);
                const auto auditEdges=assembler.sgEdgeFluxDiagnostics(state,bcs);
                (void)auditEdges;
                std::ofstream f(vela_lsa::prefix+"_direction.csv");
                f<<std::setprecision(17)<<"row,x,raw_step,capped_step,state,R,repeated_R,row_weight\n";
                for(int i=0;i<state.size();++i)f<<i<<','<<x(i)<<','<<rawStep(i)<<','<<trialStep(i)<<','<<state(i)<<','<<observed(i)<<','<<repeated(i)<<','<<activeRowWeights(i)<<'\n';
                std::ofstream m(vela_lsa::prefix+"_merit.csv");
                m<<std::setprecision(17)<<"iteration,alpha,base_norm,merit,global_e_scale,global_h_scale,scale_psi,scale_n,scale_p,weight_psi,weight_n,weight_p,residual_mode,global_mode,carrier_valid\n";
                m<<iter<<','<<alpha<<','<<residualNormFn(observed)<<','<<globalClosureLineSearchNorm(observed)<<','<<activeGlobalElectronScale<<','<<activeGlobalHoleScale<<','
                 <<residualScales.psi<<','<<residualScales.phin<<','<<residualScales.phip<<','<<residualWeights.psi<<','<<residualWeights.phin<<','<<residualWeights.phip<<','
                 <<cfg_.residualNorm<<','<<cfg_.globalContinuityClosure.mode<<','<<assembler.hasPositiveFiniteCarriers(state)<<'\n';
                if(label=="base") {
                    std::ofstream j(std::string(auditDir)+"/jacobian.csv");j<<std::setprecision(17)<<"row,column,value\n";
                    for(int k=0;k<J.outerSize();++k)for(SparseMatrixd::InnerIterator it(J,k);it;++it)j<<it.row()<<','<<it.col()<<','<<it.value()<<'\n';
                }
                vela_lsa::prefix.clear();
            };
            if(auditEnabled)auditDump(x,r,"base",0.);
'''
SEARCH_RESIDUAL=r'''
                    const VectorXd observed=assembler.residual(candidate,bcs);
                    if(auditEnabled) {
                        auditDump(candidate,observed,"trial_"+std::to_string(auditAttempt),std::ldexp(1.,-auditAttempt));
                        ++auditAttempt;
                    }
                    return observed;
'''


def build():
    a.verify(t.OUT.parent/'final_evidence.json');a.verify(t.OUT/'build_evidence.json')
    assert not (OUT/'build_evidence.json').exists()
    LOCAL.mkdir(parents=True,exist_ok=True);OUT.mkdir(parents=True,exist_ok=True)
    (LOCAL/'audit_shared.h').write_text('#pragma once\n#include <fstream>\n#include <iomanip>\n#include <string>\n#include <cstdlib>\nnamespace vela_lsa { extern std::string prefix; }\n')
    source=(t.LOCAL/'CoupledDDAssembler.cpp').read_text(encoding='utf-8')
    start=source.index('VectorXd CoupledDDAssembler::residualImpl(');end=source.index('std::vector<CoupledDDCarrierTermDiagnostic>',start)
    frag=source[start:end];assert frag.count('    return r;')==1
    source=source[:start]+frag.replace('    return r;',POISSON_DUMP+'\n    return r;')+source[end:]
    start=source.index('std::vector<CoupledDDEdgeFluxDiagnostic>\nCoupledDDAssembler::sgEdgeFluxDiagnostics(')
    end=source.index('std::vector<CoupledDDElectronTransportFactorDiagnostic>',start)
    frag=source[start:end]
    node=kernel.NODE_DUMP.replace('if(const char* auditPrefix=std::getenv("VELA_MINORITY_PRECISION_PREFIX")) {',
                                 'if(!vela_lsa::prefix.empty()) { const char* auditPrefix=vela_lsa::prefix.c_str();')
    marker='    std::vector<CoupledDDEdgeFluxDiagnostic> edges;';assert frag.count(marker)==1
    frag=frag.replace(marker,node+'\n'+marker);assert frag.count('        edges.push_back(record);')==1
    frag=frag.replace('        edges.push_back(record);',kernel.EDGE_DUMP+'\n        edges.push_back(record);')
    source=source[:start]+frag+source[end:]
    (LOCAL/'CoupledDDAssembler.cpp').write_text('#include "audit_shared.h"\n'+source,encoding='utf-8')
    # Use the very same source and compile options as the coherent core build.
    commands=a.read(c.old.LOCAL/'compile_commands.json')
    nc=list(next(x for x in commands if Path(x[x.index('-c')+1]).name=='NewtonSolver.cpp'))
    ns=Path(nc[nc.index('-c')+1]).read_text(encoding='utf-8')
    start=ns.index('        const auto runLineSearch = [&](const VectorXd& trialStep) {')
    end=ns.index('        auto ls = runLineSearch(step);',start)
    frag=ns[start:end];assert frag.count('            return lineSearch.search(')==1
    frag=frag.replace('            return lineSearch.search(',SEARCH_PREFIX+'\n            return lineSearch.search(')
    assert frag.count('                    return assembler.residual(candidate, bcs);')==1
    frag=frag.replace('                    return assembler.residual(candidate, bcs);',SEARCH_RESIDUAL)
    (LOCAL/'NewtonSolver.cpp').write_text('#include "audit_shared.h"\nnamespace vela_lsa { std::string prefix; }\n'+ns[:start]+frag+ns[end:],encoding='utf-8')
    cmds=[];replacements={}
    for name,cmd in [('CoupledDDAssembler',list(a.read(t.LOCAL/'compile_command.json'))),('NewtonSolver',nc)]:
        replacements[cmd[cmd.index('-o')+1]]=str(LOCAL/(name+'.o'))
        cmd[cmd.index('-c')+1]=str(LOCAL/(name+'.cpp'));cmd[cmd.index('-o')+1]=str(LOCAL/(name+'.o'));cmds.append(cmd)
    if (LOCAL/'compile_commands.json').exists():assert a.read(LOCAL/'compile_commands.json')==cmds
    else:a.write(LOCAL/'compile_commands.json',cmds)
    def compile_one(cmd):
        r=subprocess.run(cmd,cwd=p.REPO/'build-release',env=v.environment(),capture_output=True,encoding='utf-8',errors='replace')
        (LOCAL/(Path(cmd[cmd.index('-c')+1]).stem+'.build.log')).write_text(r.stdout+r.stderr,encoding='utf-8')
        assert r.returncode==0,r.stderr[-3500:]
    with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(compile_one,cmds))
    link=[replacements.get(x,x) for x in a.read(t.LOCAL/'link_command.json')];link[link.index('-o')+1]=str(RUNNER)
    if (LOCAL/'link_command.json').exists():assert a.read(LOCAL/'link_command.json')==link
    else:a.write(LOCAL/'link_command.json',link)
    r=subprocess.run(link,env=v.environment(),capture_output=True,encoding='utf-8',errors='replace')
    (LOCAL/'link.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr[-3000:]
    d.matrix.freeze(OUT/'build_evidence.json',[Path(__file__).resolve(),Path(kernel.__file__),t.OUT/'build_evidence.json']+[x for x in LOCAL.rglob('*') if x.is_file()]+[Path(x) for x in link if x.endswith('.o')])
    print('Built output-only final-rejection replay',flush=True)


def prepare():
    a.verify(OUT/'build_evidence.json');a.verify(prior.OUT/'validation_evidence.json')
    jobs=[];files=[Path(__file__).resolve(),OUT/'build_evidence.json',prior.OUT/'validation_evidence.json',RUNNER]
    cases={case['key']:case for case in a.read(t.OUT/'contract.json')['cases']}
    for failure in a.rows(prior.OUT/'remaining_violations.csv'):
        case=cases[failure['key']];source=prior.LOCAL/case['key']/failure['axis']/failure['label'];dest=LOCAL/case['key']
        cfg=a.read(source/'config.json');cfg['output_state_file']=str(dest/'state.csv')
        cfg['solver']['local_update_diagnostics']['csv_file']=str(dest/'updates.csv')
        assert cfg['solver'].get('damping',1.)==1.
        a.write(dest/'config.json',cfg)
        files += [dest/'config.json',source/'state.csv',source/'config.status.json',Path(cfg['state_file'])]
        jobs.append(dict(case=case,source=str(source),dest=str(dest),axis=failure['axis'],iteration=int(failure['iteration'])+1,node=int(failure['node_id'])))
    a.write(OUT/'contract.json',dict(jobs=jobs,operation='Two output-only identity replays from original restart input; dump final base, raw/capped direction, Jacobian, actual 13 backtracking candidates and physical kernel operands.',
        identity='Same output-state SHA256, accepted iteration count, failure/exit status and bit-identical repeated residual at all exported trial states.',
        precision='60/100 digit independent SG/SRH and Poisson evaluation of exact binary operands; retain rounded-product, kernel-operand, split-state and unrounded mathematical trial variants separately.',
        gates=dict(original_carrier_row=1e-6,precision_agreement_over_row_scale=1e-20,reconstruction_over_flow_scale=1e-10),
        boundary='No acceptance, algorithm, model or source-volume changes. High precision fixed-state replay is not a high precision self-consistent solve. No native DC runs.'))
    d.matrix.freeze(OUT/'freeze.json',files+[OUT/'contract.json']);print('Frozen two final-rejection identity replays',flush=True)


def run():
    a.verify(OUT/'freeze.json')
    def one(job):
        dest=Path(job['dest']);e=t.env(job['case'],job['axis'],1.)
        e.update(VELA_LS_AUDIT_DIR=str(dest),VELA_LS_AUDIT_ITERATION=str(job['iteration']))
        status=v.execute(dest/'config.json',RUNNER,e);source=Path(job['source']);old=a.read(source/'config.status.json')
        identity=a.sha(dest/'state.csv')==a.sha(source/'state.csv') and all(status[k]==old[k] for k in ('iterations','exit_code','failure_reason','converged'))
        dumps=list(dest.glob('*_direction.csv'));repeat=all(all(row['R']==row['repeated_R'] for row in a.rows(f)) for f in dumps)
        row=dict(key=job['case']['key'],identity=identity,repeat_residual_identity=repeat,snapshots=len(dumps),iteration=status['iterations'],failure=status['failure_reason'])
        print(row,flush=True);assert identity and repeat and len(dumps)==14,row
        return row
    with ThreadPoolExecutor(max_workers=2) as pool:rows=list(pool.map(one,a.read(OUT/'contract.json')['jobs']))
    a.write_csv(OUT/'identity.csv',rows)
    d.matrix.freeze(OUT/'replay_evidence.json',[OUT/'freeze.json',OUT/'identity.csv']+[x for j in a.read(OUT/'contract.json')['jobs'] for x in Path(j['dest']).rglob('*') if x.is_file()])


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('build','prepare','run'));globals()[parser.parse_args().action]()
