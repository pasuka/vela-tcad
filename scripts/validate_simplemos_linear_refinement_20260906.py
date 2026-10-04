"""Bounded self-consistent A/B: four linear residual corrections, unchanged Newton gates."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import copy
import json
from pathlib import Path
import subprocess
import time
import math
import audit_simplemos_minority_residual_20260906 as v

a=v.a;d=v.d
LOCAL=d.REPO/'build-release/simplemos_linear_refinement_20260906'
OUT=d.ROOT/'linear_refinement_20260906'
RUNNER=LOCAL/'runner.exe'
HEADER=d.REPO/'scripts/diagnostics/simplemos_linear_refinement.hpp'
TEST=d.REPO/'tests/diagnostics/test_simplemos_linear_refinement.cpp'

PATCH=r'''
            if (const char* prefix=std::getenv("VELA_SIMPLEMOS_LINEAR_REFINEMENT_PREFIX")) {
                const auto terms=assembler.carrierContinuityEquationTermDiagnostics(x,bcs);
                VectorXd scales=VectorXd::Zero(3*N);int active=0,zero=0;
                for(int b=1;b<=2;++b)for(int i=0;i<N;++i) {
                    const auto& t=terms[i];if(!(b==1?t.electronContinuityActive:t.holeContinuityActive))continue;
                    ++active;
                    const double s=std::max({b==1?t.electronFluxAbsSum:t.holeFluxAbsSum,
                        std::abs(b==1?t.electronRecombination:t.holeRecombination),
                        std::abs(b==1?t.electronImpact:t.holeImpact)});
                    scales(b*N+i)=s;if(s==0)++zero;
                }
                std::ofstream summary(std::string(prefix)+"_summary.csv",iter==1?std::ios::out:std::ios::app);
                std::ofstream selected(std::string(prefix)+"_selected.csv",iter==1?std::ios::out:std::ios::app);
                if(!summary||!selected)throw std::runtime_error("Cannot write isolated refinement diagnostics");
                summary<<std::setprecision(17);selected<<std::setprecision(17);
                if(iter==1) {
                    summary<<"iteration,correction,active_rows,zero_scale_rows,max_carrier_linear_ratio,rows_above_1e_minus_8,worst_row,componentwise_backward_max,physical_defect_l2,physical_rhs_l2,step_inf_V\n";
                    selected<<"iteration,correction,carrier,node,row_scale,original_residual,linear_defect_hp,step_V\n";
                }
                const VectorXd weights=cfg_.continuityRowScaling.enabled?activeRowWeights:VectorXd::Ones(3*N);
                step=simplemos_diagnostic::refine(J,r,weights,step,linearSolver,4,
                    [&](int k,const VectorXd& direction,const simplemos_diagnostic::Defect& defect) {
                        double maximum=0;int count=0,worst=-1;
                        for(int row=N;row<3*N;++row)if(scales(row)>0) {
                            const double ratio=std::abs(defect.values(row))/scales(row);
                            if(ratio>maximum){maximum=ratio;worst=row;}
                            if(ratio>1e-8)++count;
                        }
                        summary<<iter<<','<<k<<','<<active<<','<<zero<<','<<maximum<<','<<count<<','<<worst<<','<<defect.backwardMax<<','<<defect.values.norm()<<','<<r.norm()<<','<<direction.lpNorm<Eigen::Infinity>()*potentialScale<<'\n';
                        for(int b=1;b<=2;++b)for(int node:{41,453,973,949,355})if(node<N) {
                            int row=b*N+node;
                            selected<<iter<<','<<k<<','<<(b==1?"electron":"hole")<<','<<node<<','<<scales(row)<<','<<r(row)<<','<<defect.values(row)<<','<<direction(row)*potentialScale<<'\n';
                        }
                    });
            }
'''


def env(refined=False,path=None):
    e=v.m.env()
    for key in ('VELA_MINORITY_PRECISION_PREFIX','VELA_MINORITY_TRIAL_CSV','VELA_SIMPLEMOS_LINEAR_REFINEMENT_PREFIX'):
        e.pop(key,None)
    if refined:e['VELA_SIMPLEMOS_LINEAR_REFINEMENT_PREFIX']=str(path/'linear')
    return e


def build():
    a.verify(v.OUT/'validation_evidence.json');LOCAL.mkdir(parents=True,exist_ok=False)
    source=(v.LOCAL/'NewtonSolver.cpp').read_text()
    marker='''            } else {
                step = linearSolver.solve(J, -r);
            }
        } catch (const std::runtime_error&) {'''
    assert source.count(marker)==1
    source=source.replace(marker,marker.replace('        } catch',PATCH+'        } catch'))
    (LOCAL/'NewtonSolver.cpp').write_text('#include "simplemos_linear_refinement.hpp"\n'+source,newline='\n')
    original=a.read(v.LOCAL/'build_command.json')
    args=[str(LOCAL/'NewtonSolver.cpp') if p==str(v.LOCAL/'NewtonSolver.cpp') else p for p in original]
    args.insert(1,'-I'+str(HEADER.parent));args[-1]=str(RUNNER);a.write(LOCAL/'build_command.json',args)
    testargs=[p for p in args if not p.endswith('.cpp')]
    testargs[-1]=str(LOCAL/'test_refinement.exe')
    index=next(i for i,p in enumerate(testargs) if p.endswith('libvela_core.a'))
    testargs.insert(index,str(TEST));index=testargs.index('-o')
    testargs[index:index]=['D:/msys64/ucrt64/lib/libCatch2Main.a','D:/msys64/ucrt64/lib/libCatch2.a']
    a.write(LOCAL/'test_build_command.json',testargs)
    for label,command in (('test_build',testargs),('build',args)):
        p=subprocess.run(command,env=env(),capture_output=True,text=True)
        (LOCAL/(label+'.log')).write_text(p.stdout+p.stderr);assert p.returncode==0,p.stderr[-4000:]
        print(label,'passed',flush=True)
        if label=='test_build':
            t=subprocess.run([str(LOCAL/'test_refinement.exe')],env=env(),capture_output=True,text=True)
            (LOCAL/'test_result.log').write_text(t.stdout+t.stderr);assert t.returncode==0,t.stdout+t.stderr
            print(t.stdout,flush=True)


def prepare():
    a.verify(v.OUT/'validation_evidence.json')
    def write_config(path,cfg):
        if path.exists():assert a.read(path)==cfg,path
        else:a.write(path,cfg)
    files=[Path(__file__).resolve(),HEADER,TEST,RUNNER,LOCAL/'test_refinement.exe',LOCAL/'NewtonSolver.cpp',
        LOCAL/'build_command.json',LOCAL/'test_result.log',v.OUT/'validation_evidence.json']
    files += [Path(p) for p in a.read(LOCAL/'build_command.json') if p.endswith(('.cpp','.a'))]
    jobs=[];baselines=[]
    previous=a.read(v.m.OUT/'contract.json')['jobs']
    for trace in a.read(v.OUT/'contract.json')['traces']:
        src=Path(trace['source']);job=next(j for j in previous if Path(j['config']).parent==src)
        # Full original initializations, including hole_plus: no saved-state one-step restart here.
        for mode in ('refined','disabled'):
            sentinel=job['initialization']=='vela' and job['vd']==.05
            if mode=='disabled' and not sentinel:continue
            dest=LOCAL/mode/Path(trace['tag']);cfg=a.read(src/'config.json');cfg['output_state_file']=str(dest/'state.csv')
            cfg['solver']['local_update_diagnostics']=dict(enabled=True,nodes=list(range(1480 if job['device']=='n19' else 1482)),
                csv_file=str(dest/'updates.csv'),first_iterations=200,every_iterations=1)
            assert cfg['solver']['max_iter']==cfg['solver']['carrier_row_convergence']['min_newton_max_iter']==200
            assert cfg['solver']['carrier_row_convergence']['eps_row']==1e-6
            assert cfg['solver'].get('carrier_regularization_scale',0)==0
            check=copy.deepcopy(cfg);check['output_state_file']=a.read(src/'config.json')['output_state_file']
            check['solver'].pop('local_update_diagnostics');assert check==a.read(src/'config.json')
            write_config(dest/'config.json',cfg)
            entry={**job,'tag':trace['tag'],'source':str(src),'config':str(dest/'config.json'),'mode':mode}
            (jobs if mode=='refined' else baselines).append(entry)
            files += [dest/'config.json',src/'config.json',src/'state.csv',src/'config.status.json',src/'result.json',Path(cfg['state_file'])]
            for label,base in (('all_row',cfg),('legacy',a.read(Path(job['legacy_config'])))):
                probe=copy.deepcopy(base);probe.pop('output_state_file',None);probe['solver'].pop('local_update_diagnostics',None)
                probe.update(simulation_type='newton_carrier_term_probe',state_file=str(dest/'state.csv'),output_csv=str(dest/(label+'.csv')),
                    carrier_term_probe={'solved_equation_terms':True})
                probe['solver']['carrier_row_convergence']['mode']='report'
                probe['solver']['global_continuity_closure']=dict(mode='enforce',tolerance=1e-6,source_floor=1e-10)
                write_config(dest/(label+'.json'),probe);files.append(dest/(label+'.json'))
    assert len(jobs)==7 and len(baselines)==2
    write_config(OUT/'contract.json',dict(status='frozen_before_execution',jobs=jobs,disabled_sentinels=baselines,
        single_change='After the original SparseLU raw solve, apply exactly four same-matrix corrections with 100-digit accumulated Jdx+F and the SAME row weights. No equilibration, no change to caps, line-search, physical equations, Jacobian or final gates.',
        linear_monitor='At every Newton iteration check all active continuity rows using the current pre-step physical flux/source scale; 1e-8 is a diagnostic target, 100x below the unchanged 1e-6 nonlinear gate. Also report componentwise backward error over all matrix rows. Neither monitor changes acceptance.',
        comparison='Seven frozen original initializations vs their retained baseline outputs; two disabled-binary sentinels must reproduce original state hashes/iterations/exit/failure before treatment runs.',
        qualification='Original nonlinear exit AND all 1814 active rows with no zero-scale exclusion AND original global closure AND KCL/Id<=1e-8. Failed states remain failed regardless of current agreement.',
        four_seed_invariance='Only the selected n23 Vd=.05 Vg=1 V case has all four initializations. Apply the prior frozen phi 1e-6 V, density 1e-4, Id 1e-6 invariance criteria only if all four qualify.',
        backend='Eigen SparseLU UCRT64 Release',environment={k:x for k,x in env().items() if k.startswith('VELA_')},
        remote_runs=0,production_changes=False,m82_released=False,m83_released=False))
    files.append(OUT/'contract.json')
    a.write(OUT/'freeze.json',{'input_hashes':{a.rel(p) if p.is_relative_to(d.REPO) else p.as_posix():a.sha(p) for p in sorted(set(files))}})
    print('Frozen 7 refinement solves and 2 disabled sentinels',flush=True)


def execute(path,refined=False):
    target=path.with_suffix('.status.json')
    if target.exists():return a.read(target)
    start=time.monotonic();p=subprocess.run([str(RUNNER),'--config',str(path),'--log','off'],env=env(refined,path.parent),capture_output=True,text=True)
    path.with_suffix('.stdout.txt').write_text(p.stdout);path.with_suffix('.stderr.txt').write_text(p.stderr)
    s=json.loads(p.stdout.strip().splitlines()[-1]);s.update(exit_code=p.returncode,elapsed_seconds=time.monotonic()-start)
    a.write(target,s);return s


def one(job):
    path=Path(job['config']);dest=path.parent;s=execute(path,job['mode']=='refined')
    assert (dest/'state.csv').exists(),s
    if job['mode']=='disabled':
        old=a.read(Path(job['source'])/'config.status.json')
        assert a.sha(dest/'state.csv')==a.sha(Path(job['source'])/'state.csv')
        assert all(s[k]==old[k] for k in ('iterations','converged','exit_code','failure_reason'))
    audits={label:execute(dest/(label+'.json')) for label in ('all_row','legacy')}
    current=s['contact_currents_A_per_um']['drain'];kcl=abs(math.fsum(s['contact_currents_A_per_um'].values()))/max(abs(current),1e-300)
    result=dict(**job,iterations=s['iterations'],converged=s['converged'],exit_code=s['exit_code'],failure_reason=s['failure_reason'],
        current_A_per_um=current,kcl_over_Id=kcl,elapsed_seconds=s['elapsed_seconds'])
    for label,probe in audits.items():
        gate=probe['carrier_row_convergence'];assert gate['qualified_row_count']==1814 if label=='all_row' else True
        result.update({label+'_violations':gate['violation_count'],label+'_max_ratio':gate['max_ratio'],
            label+'_rows_satisfied':gate['satisfied'],label+'_global_satisfied':probe['global_continuity_closure']['satisfied'],
            label+'_qualified':s['converged'] and s['exit_code']==0 and probe['exit_code']==0 and gate['satisfied'] and probe['global_continuity_closure']['satisfied'] and kcl<=1e-8})
    a.write(dest/'result.json',result)
    print(job['mode'],job['tag'],s['iterations'],s['failure_reason'],'qualified',result['all_row_qualified'],'rows',result['all_row_violations'],flush=True)
    return result


def run():
    a.verify(OUT/'freeze.json');c=a.read(OUT/'contract.json')
    baseline=[one(j) for j in c['disabled_sentinels']]
    a.write_csv(OUT/'disabled_sentinels.csv',baseline)
    with ThreadPoolExecutor(max_workers=2) as pool:rows=list(pool.map(one,c['jobs']))
    a.write_csv(OUT/'runs.csv',rows)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('build','prepare','run'));globals()[p.parse_args().action]()
