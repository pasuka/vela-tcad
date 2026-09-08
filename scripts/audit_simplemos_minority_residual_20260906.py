"""Frozen output-only precision probes and selected Newton rejection replays."""
import argparse
import copy
import json
from pathlib import Path
import subprocess
import validate_simplemos_minority_invariance_20260906 as m
import trace_simplemos_strict_rejection as oldtrace

a=m.a;d=m.d
LOCAL=d.REPO/'build-release/simplemos_minority_residual_20260906'
OUT=d.ROOT/'minority_residual_20260906';RUNNER=LOCAL/'runner.exe'

NODE_DUMP=r'''
    std::ofstream auditEdges;
    if(const char* auditPrefix=std::getenv("VELA_MINORITY_PRECISION_PREFIX")) {
        auditEdges.open(std::string(auditPrefix)+"_edges.csv");
        auditEdges<<std::setprecision(17)<<"edge,node0,node1,ni0,ni1,epsi0,epsi1,ephi0,ephi1,hpsi0,hpsi1,hphi0,hphi1,ncoef,pcoef,scale,nflux,pflux,Vt\n";
        std::ofstream nodes(std::string(auditPrefix)+"_nodes.csv");
        nodes<<std::setprecision(17)<<"node,psi,eref,einc,href,hinc,n,p,nsrh,ni,dphi,taun,taup,volume,source_factor,scale,rate,Vt\n";
        const Real sf=scaling_.enabled?scaling_.unitSystem.continuitySourceIntegralFactor():1.0;
        for(Index i=0;i<Nidx;++i) {
            Real ei=x(phinOffset()+i)*potentialScale,hi=x(phipOffset()+i)*potentialScale;
            Real er=electronQuasiFermiReferenceAt(i),hr=holeQuasiFermiReferenceAt(i);
            Real dp=referencedDifference(er,ei,hr,hi);
            Real ns=electronSrhDensityAt(i,psi(i)-er,ei);
            Real concentration=recombination_.srhDopingConcentration(doping_.donors(i),doping_.acceptors(i));
            Real rate=ni_[i]>0?nodeRecombinationRate(i,n(i),ns,p(i),dp):0.0;
            nodes<<i<<','<<psi(i)<<','<<er<<','<<ei<<','<<hr<<','<<hi<<','<<n(i)<<','<<p(i)<<','<<ns<<','<<ni_[i]<<','<<dp<<','
                 <<recombination_.electronLifetime(concentration)<<','<<recombination_.holeLifetime(concentration)<<','<<vol_[i]<<','<<sf<<','<<continuityScale<<','<<rate<<','<<Vt_<<'\n';
        }
    }
'''
EDGE_DUMP=r'''
        if(auditEdges.is_open()) {
            auditEdges<<e<<','<<i<<','<<j<<','<<ni_[idxI]<<','<<ni_[idxJ]<<','<<electronPsiRelative_i<<','<<electronPsiRelative_j<<','<<phin_i<<','<<phin_j<<','
                <<holePsiRelative_i<<','<<holePsiRelative_j<<','<<phip_i<<','<<phip_j<<','
                <<mun*Vt_*fieldFactor*couple_[e]/h<<','<<mup*Vt_*fieldFactor*couple_[e]/h<<','<<continuityScale<<','<<record.electronFlux<<','<<record.holeFlux<<','<<Vt_<<'\n';
        }
'''
TRIAL_DUMP=r'''
                    const VectorXd evaluated=assembler.residual(candidate,bcs);
                    if(const char* path=std::getenv("VELA_MINORITY_TRIAL_CSV")) {
                        const bool header=!std::filesystem::exists(path);
                        std::ofstream out(path,std::ios::app);out<<std::setprecision(17);
                        if(header)out<<"iteration,attempt,node,psi,eref,einc,href,hinc,rpsi,rn,rp\n";
                        for(int i=0;i<N;++i)out<<iter<<','<<auditTrial<<','<<i<<','<<candidate(i)*potentialScale<<','
                            <<assembler.electronQuasiFermiReferenceAt(i)<<','<<candidate(N+i)*potentialScale<<','
                            <<assembler.holeQuasiFermiReferenceAt(i)<<','<<candidate(2*N+i)*potentialScale<<','
                            <<evaluated(i)<<','<<evaluated(N+i)<<','<<evaluated(2*N+i)<<'\n';
                    }
                    ++auditTrial;
                    return evaluated;
'''


def build():
    a.verify(m.OUT/'validation_evidence.json');LOCAL.mkdir(parents=True,exist_ok=False)
    src=(m.surface.LOCAL/'CoupledDDAssembler.cpp').read_text()
    start=src.index('std::vector<CoupledDDEdgeFluxDiagnostic>\nCoupledDDAssembler::sgEdgeFluxDiagnostics(')
    end=src.index('std::vector<CoupledDDElectronTransportFactorDiagnostic>',start)
    frag=src[start:end];marker='    std::vector<CoupledDDEdgeFluxDiagnostic> edges;';assert frag.count(marker)==1
    frag=frag.replace(marker,NODE_DUMP+'\n'+marker)
    marker='        edges.push_back(record);';assert frag.count(marker)==1;frag=frag.replace(marker,EDGE_DUMP+marker)
    (LOCAL/'CoupledDDAssembler.cpp').write_text('#include <fstream>\n'+src[:start]+frag+src[end:],newline='\n')
    ns=(d.REPO/'build-release/simplemos_convergence_audit/NewtonSolver.cpp').read_text()
    start=ns.index('        const auto runLineSearch = [&](const VectorXd& trialStep) {')
    end=ns.index('        auto ls = runLineSearch(step);',start)
    frag=ns[start:end].replace('return lineSearch.search(','int auditTrial=0;\n            return lineSearch.search(',1)
    marker='                    return assembler.residual(candidate, bcs);';assert frag.count(marker)==1
    ns=ns[:start]+frag.replace(marker,TRIAL_DUMP)+ns[end:]
    (LOCAL/'NewtonSolver.cpp').write_text('#include <cstdlib>\n'+ns,newline='\n')
    runner=(m.surface.m.LOCAL/'runner.cpp').read_text();start=runner.index('nlohmann::json runNewtonSolveFromState(');end=runner.index('std::vector<vela::Real> readNodeScalarCsv',start)
    frag=runner[start:end];assert frag.count('    return {')==1
    frag=frag.replace('    return {',oldtrace.TRACE+'\n    return {\n        {"newton_trace", trace},',1)
    (LOCAL/'runner.cpp').write_text(runner[:start]+frag+runner[end:],newline='\n')
    args=[str(LOCAL/Path(s).name) if Path(s).name in ('runner.cpp','CoupledDDAssembler.cpp','NewtonSolver.cpp') else s for s in a.read(m.surface.LOCAL/'build_command.json')]
    args[-1]=str(RUNNER);a.write(LOCAL/'build_command.json',args)
    p=subprocess.run(args,env=m.env(),capture_output=True,text=True);(LOCAL/'build.log').write_text(p.stdout+p.stderr)
    assert p.returncode==0,p.stderr[-3000:];print('Built output-only precision and trace overlay',flush=True)


def prepare():
    a.verify(m.OUT/'validation_evidence.json');files=[Path(__file__).resolve(),m.OUT/'validation_evidence.json',RUNNER]
    files += [d.REPO/'scripts/analyze_simplemos_minority_residual_20260906.py',d.REPO/'tests/regression/test_simplemos_minority_residual_precision.py']
    files += [LOCAL/n for n in ('runner.cpp','CoupledDDAssembler.cpp','NewtonSolver.cpp','build_command.json')]
    probes=[];traces=[]
    for row in a.rows(m.OUT/'qualification.csv'):
        original=next(j for j in a.read(m.OUT/'contract.json')['jobs'] if j['case']==row['case'] and float(j['vg'])==float(row['vg']) and j['initialization']==row['initialization'])
        source=Path(original['config']).parent;tag=source.relative_to(m.LOCAL);dest=LOCAL/'probes'/tag
        cfg=a.read(source/'config.json');cfg.pop('output_state_file');cfg.update(simulation_type='sg_edge_flux_probe',state_file=str(source/'state.csv'),output_csv=str(dest/'edges.csv'))
        a.write(dest/'config.json',cfg);probes.append({**original,'config':str(dest/'config.json'),'source':str(source),'tag':str(tag)})
        files += [dest/'config.json',source/'state.csv',source/'all_row.csv',source/'all_row.status.json']
        # Identity replays: both n19 failing seeds; both n23 passing seeds; a first-trial rejection; a high-Vd control.
        replay=(row['device']=='n19' and float(row['vd'])==.05 and float(row['vg'])==.35 and row['initialization'] in ('vela','native')) or (row['device']=='n23' and float(row['vd'])==.05 and float(row['vg'])==1 and row['initialization'] in ('vela','native','hole_minus')) or (row['device']=='n23' and float(row['vd'])==1 and float(row['vg'])==.35 and row['initialization']=='vela')
        single=row['device']=='n23' and float(row['vd'])==.05 and float(row['vg'])==1 and row['initialization']=='hole_plus'
        if not (replay or single):continue
        dest=LOCAL/'traces'/tag;cfg=a.read(source/'config.json');cfg['output_state_file']=str(dest/'state.csv');nodes=list(range(1480 if row['device']=='n19' else 1482))
        cfg['solver']['local_update_diagnostics']=dict(enabled=True,nodes=nodes,csv_file=str(dest/'updates.csv'),first_iterations=200,every_iterations=1)
        if single:
            cfg['state_file']=str(source/'state.csv');cfg['solver']['max_iter']=1;cfg['solver']['carrier_row_convergence']['min_newton_max_iter']=1
        a.write(dest/'config.json',cfg);files.append(dest/'config.json');traces.append(dict(tag=str(tag),source=str(source),config=str(dest/'config.json'),mode='one_step_saved_failure' if single else 'identity_replay'))
    assert len(probes)==20 and len(traces)==7
    a.write(OUT/'contract.json',dict(status='frozen_before_execution',probes=probes,traces=traces,
        purpose='Separate kernel arithmetic, QF-coordinate conversion and nonclosed state; inspect raw/capped Newton updates and actual line search.',
        precision='Independent Decimal 60/100 digit references, exact binary-float inputs. Mobility and lifetimes fixed at production values. Kernel-rounded and split-coordinate references kept separately.',
        gates=dict(row_eps=1e-6,reconstruction_scale_relative=1e-10,precision_agreement_scale_relative=1e-20),
        identity='Six output-only original-initialization replays must match final-state SHA256, iteration count and exit status. One saved max-iteration failure gets only one diagnostic step; it is not a convergence qualification or reconstruction of the historical 200th step.',
        unchanged='Physical equations, original tolerances, references, production sources and earlier failure records unchanged. No new Sentaurus runs.'))
    files.append(OUT/'contract.json');d.matrix.freeze(OUT/'freeze.json',files)
    print('Frozen 20 probes, 6 identity replays and 1 single-step diagnostic',flush=True)


def execute(path,mode):
    status=path.with_suffix('.status.json')
    if status.exists():return a.read(status)
    e=m.env()
    if mode=='probe':e['VELA_MINORITY_PRECISION_PREFIX']=str(path.parent/'kernel')
    else:e['VELA_MINORITY_TRIAL_CSV']=str(path.parent/'trials.csv')
    p=subprocess.run([str(RUNNER),'--config',str(path),'--log','off'],env=e,capture_output=True,text=True)
    path.with_suffix('.stdout.txt').write_text(p.stdout);path.with_suffix('.stderr.txt').write_text(p.stderr)
    s=json.loads(p.stdout.strip().splitlines()[-1]);s['exit_code']=p.returncode;a.write(status,s);return s


def probes():
    a.verify(OUT/'freeze.json')
    for job in a.read(OUT/'contract.json')['probes']:
        s=execute(Path(job['config']),'probe');assert s['exit_code']==0,s
        print(job['tag'],'probe complete',flush=True)


def traces():
    a.verify(OUT/'freeze.json');results=[]
    for job in a.read(OUT/'contract.json')['traces']:
        path=Path(job['config']);s=execute(path,'trace');source=Path(job['source']);old=a.read(source/'config.status.json')
        identical=a.sha(path.parent/'state.csv')==a.sha(source/'state.csv')
        if job['mode']=='identity_replay':assert identical and s['iterations']==old['iterations'] and s['exit_code']==old['exit_code'] and s['failure_reason']==old['failure_reason'],job
        results.append(dict(**job,state_identical=identical,iterations=s['iterations'],failure=s['failure_reason'],converged=s['converged']))
        print(job['tag'],job['mode'],s['iterations'],s['failure_reason'],'identical',identical,flush=True)
    a.write_csv(OUT/'traces.csv',results)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('build','prepare','probes','traces'));globals()[p.parse_args().action]()
