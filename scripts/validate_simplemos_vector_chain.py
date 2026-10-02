"""Isolated vector-drive Jacobian chain completion and fixed-state validation."""
import argparse
import copy
import json
import math
from pathlib import Path
import subprocess
import numpy as np
import check_simplemos_actual_direction_jvp as jvp
import seal_simplemos_actual_direction_jvp as previous_seal

p=jvp.prior
REPO,ROOT=p.REPO,p.ROOT
LOCAL=REPO/'build-release/simplemos_vector_chain_validation'
OUT=ROOT/'vector_chain_validation'
RUNNER=LOCAL/'chain_runner.exe'
FREEZE=OUT/'freeze.json'
CONTRACT=OUT/'contract.json'
FLAG='VELA_VALIDATE_VECTOR_CHAIN_ENABLE'

CHAIN=r'''
        // Isolated diagnostic: add only the missing field-drive chain term.
        if (std::getenv("VELA_VALIDATE_VECTOR_CHAIN_ENABLE") && vectorQfMobility &&
            mobilityConfig_.jacobianFieldDerivatives) {
            for (int carrier=0;carrier<2;++carrier) {
                const bool electron=carrier==0;
                const Real baseMobility=electron?mun:mup;
                if (!(baseMobility>0.0)) continue;
                const VectorXd& qf=electron?phinState:phipState;
                const Real baseFlux=electron
                    ? edgeElectronTransportFlux(e,i,j,h,psi_i,psi_j,phin_i,phin_j,phip_i,phip_j,-1.0,mesh_.numNodes(),0.0)
                    : edgeHoleTransportFlux(e,i,j,h,psi_i,psi_j,phip_i,phip_j,phin_i,phin_j,-1.0,mesh_.numNodes(),0.0);
                const int offset=electron?phinOffset():phipOffset();
                for(std::size_t k=0;k<edgeKernel.avalancheStencilNodeCount;++k) {
                    const Index node=edgeKernel.avalancheStencilNodes[k];
                    const Real step=1e-7;
                    const Real ep=transportVectorMobilityField(e,[&](Index t) {return qf(t)+(t==node?step:0.0);});
                    const Real em=transportVectorMobilityField(e,[&](Index t) {return qf(t)-(t==node?step:0.0);});
                    const auto type=electron?CarrierType::Electron:CarrierType::Hole;
                    const Real mp=cachedEdgeMobility(e,type,ep,&psi,n(i),n(j),p(i),p(j));
                    const Real mm=cachedEdgeMobility(e,type,em,&psi,n(i),n(j),p(i),p(j));
                    const Real derivative=baseFlux/baseMobility*(mp-mm)/(2*step);
                    add(offset+i,offset+static_cast<int>(node),derivative);
                    add(offset+j,offset+static_cast<int>(node),-derivative);
                }
            }
        }
'''


def env(enabled=True,charge=None):
    e=p.env(charge);e.pop(FLAG,None)
    if enabled:e[FLAG]='1'
    return e


def build():
    if RUNNER.exists():raise FileExistsError(RUNNER)
    LOCAL.mkdir(parents=True,exist_ok=True)
    source=(REPO/'src/equation/CoupledDDAssembler.cpp').read_text()
    marker='        if (sgCurrentAvalanche && !cellLocalAvalanche &&\n            impactIonizationConfig_.sourceJacobianMode == "finite_difference") {'
    assert source.count(marker)==1
    source=source.replace(marker,CHAIN+'\n'+marker)
    marker='        recombination_.bandToBandEnabled() || surfaceMobilityEnabled_;'
    assert source.count(marker)==1
    source=source.replace(marker,'        recombination_.bandToBandEnabled() || surfaceMobilityEnabled_ ||\n        (std::getenv("VELA_VALIDATE_VECTOR_CHAIN_ENABLE") && vectorQfMobility);')
    marker='    carrierStatisticsModel_ = carrierStatisticsModel(carrierStatistics_);'
    assert source.count(marker)==1
    source=source.replace(marker,p.HOOK+'\n'+marker)
    (LOCAL/'CoupledDDAssembler.cpp').write_text('#include <cstdlib>\n#include <fstream>\n'+source,newline='\n')
    text=(REPO/'build-release/m79b_numerical_calibration/vela_example_runner_calibrated.cpp').read_text()
    text=text.replace('int main(int argc, char** argv)',jvp.FUNCTION+'\nint main(int argc, char** argv)')
    marker='        } else if (type == "newton_residual_probe") {'
    assert text.count(marker)==1
    text=text.replace(marker,'        } else if (type == "actual_direction_jvp") {\n            status.update(runActualDirection(configFile, cfg));\n'+marker)
    (LOCAL/'runner.cpp').write_text(text,newline='\n')
    args=p.audit.read(REPO/'build-release/simplemos_strict_rejection_trace/build_command.json')
    args=[str(LOCAL/'runner.cpp') if a.endswith('runner.cpp') else a for a in args]
    args.insert(next(i for i,a in enumerate(args) if a.endswith('libvela_core.a')),str(LOCAL/'CoupledDDAssembler.cpp'))
    args[-1]=str(RUNNER)
    p.audit.write(LOCAL/'build_command.json',args)
    proc=subprocess.run(args,env=env(False),capture_output=True,text=True)
    (LOCAL/'build.log').write_text(proc.stdout+proc.stderr)
    if proc.returncode:raise RuntimeError(proc.stderr[-2000:])
    print('Built isolated chain correction; production files untouched',flush=True)


def prepare():
    previous_seal.verify()
    paths=[Path(__file__),p.FREEZE,jvp.FREEZE,jvp.OUT/'evidence.json',p.OUT/'evidence.json',
           p.CONTRACT,p.OUT/'calibration.csv',p.BASE/'state.csv',REPO/'src/equation/CoupledDDAssembler.cpp']
    for mode in ('off','on'):
        for case in ('1e13','5e12'):
            cfg=p.audit.read(jvp.LOCAL/case/'config.json')
            cfg['output_csv']=str(LOCAL/mode/case/'jvp.csv')
            dest=LOCAL/mode/case/'config.json';p.audit.write(dest,cfg)
            paths += [dest]+[Path(v) for v in cfg['direction_files'].values()]
            paths += [Path(cfg[k]) for k in ('mesh_file','materials_file','node_doping_file','state_file','plus_state_file','minus_state_file')]
    cfg=p.audit.read(p.native.upstream.LOCAL/'m65_n23_vd_0p050000_endpoint/config.json')
    cfg['output_csv']=str(LOCAL/'adjoint/adjoint.csv')
    cfg['electron_volume_response']['output_csv']=str(LOCAL/'adjoint/response.csv')
    p.audit.write(LOCAL/'adjoint/config.json',cfg)
    for name,_ in p.CASES:
        cfg=p.audit.read(p.LOCAL/name/'config.json')
        cfg['output_state_file']=str(LOCAL/'dc'/name/'state.csv')
        p.audit.write(LOCAL/'dc'/name/'config.json',cfg)
        paths.append(p.LOCAL/name/'charge.csv')
    paths += [x for x in LOCAL.rglob('*') if x.is_file()]
    for a in p.audit.read(LOCAL/'build_command.json'):
        if Path(a).is_file():paths.append(Path(a))
    p.audit.write(CONTRACT,{'status':'frozen_before_execution','single_axis':'Add F/mu * dmu/d(vector_qf_field) * d(vector_qf_field)/d(qf_node) in existing full adjacent-cell stencil, for both carriers; equal/opposite endpoint rows; existing state-density and surface derivatives retained.',
        'field_chain_finite_difference_step_V':1e-7,'switch':FLAG,'gates':{'residual_rows_identical':True,'old_Jv_rows_identical_when_off':True,'weighted_defect_reduction_minimum':.99,'fd_adjoint_relative':.001,'base_current_relative':1e-8,'strict_KCL_over_Id':1e-8},
        'release_order':'Fixed-state Jv and current invariance, corrected adjoint versus already computed nonlinear FD, then five strict reclosures if both pass.',
        'scope':'n23 low Vd only. No production default, current functional, physical residual or Sentaurus change. No M82/M83 release.'})
    paths.append(CONTRACT)
    p.audit.write(FREEZE,{'input_hashes':{(p.audit.rel(x) if x.is_relative_to(REPO) else str(x)):p.audit.sha(x) for x in sorted(set(paths))}})


def execute(cfgpath,enabled=True,charge=None):
    proc=subprocess.run([str(RUNNER),'--config',str(cfgpath),'--log','off'],env=env(enabled,charge),capture_output=True,text=True)
    cfgpath.with_suffix('.stdout.txt').write_text(proc.stdout);cfgpath.with_suffix('.stderr.txt').write_text(proc.stderr)
    status=json.loads(proc.stdout.strip().splitlines()[-1]);status['exit_code']=proc.returncode
    p.audit.write(cfgpath.with_suffix('.status.json'),status)
    if proc.returncode:raise RuntimeError(str(status))
    return status


def check():
    p.audit.verify(FREEZE)
    adj={int(r['node_id']):r for r in p.audit.rows(p.native.upstream.LOCAL/'m65_n23_vd_0p050000_endpoint/adjoint.csv')}
    w={'psi':'lambda_poisson','phin':'lambda_electron','phip':'lambda_hole'}
    checks=[]
    for case in ('1e13','5e12'):
        before=p.audit.rows(jvp.LOCAL/case/'jvp.csv')
        for mode in ('off','on'):
            execute(LOCAL/mode/case/'config.json',mode=='on')
            after=p.audit.rows(LOCAL/mode/case/'jvp.csv')
            assert len(after)==len(before)
            sameF=all(a['fd']==b['fd'] and a['actual_endpoint_fd']==b['actual_endpoint_fd'] for a,b in zip(after,before))
            sameJ=all(a['analytic']==b['analytic'] for a,b in zip(after,before))
            assert sameF and (mode=='on' or sameJ)
            for h in (1,.5,.25):
                group=[(a,b) for a,b in zip(after,before) if a['mode']=='all' and float(a['step'])==h]
                sums=[]
                for pos in (0,1):
                    sums.append(math.fsum(float(adj[int(pair[pos]['node_id'])][w[pair[pos]['row_block']]])*(float(pair[pos]['analytic'])-float(pair[pos]['fd'])) for pair in group))
                reduction=1-abs(sums[0]/sums[1])
                checks.append({'case':case,'mode':mode,'step':h,'residual_identical':sameF,'Jv_identical':sameJ,'old_weighted_defect_A_per_um':sums[1],'new_weighted_defect_A_per_um':sums[0],'defect_reduction':reduction,'pass':sameF and (reduction>=.99 if mode=='on' else sameJ)})
    p.audit.write_csv(OUT/'jvp_checks.csv',checks)
    assert all(r['pass'] for r in checks)
    status=execute(LOCAL/'adjoint/config.json')
    newadj=p.audit.rows(LOCAL/'adjoint/adjoint.csv')
    charge=p.audit.rows(p.LOCAL/'plus_1e13/charge.csv')
    factor=p.audit.read(p.CONTRACT)['frozen_poisson_residual_per_physical_charge']
    pred=math.fsum(float(a['lambda_poisson'])*factor*float(q['charge_C_per_m']) for a,q in zip(newadj,charge))
    actual=float(p.audit.rows(p.OUT/'calibration.csv')[0]['central_delta_A_per_um'])
    base=float(p.audit.rows(p.OUT/'case_ledger.csv')[0]['current_A_per_um'])
    result={'Jv_checks_passed':True,'prediction_A_per_um_at_1e13':pred,'old_self_consistent_FD_A_per_um':actual,'relative_error':abs(pred/actual-1),
        'base_current_relative_change':abs(status['current_A_per_um']/base-1),'adjoint_relative_residual':status['adjoint_relative_residual'],
        'minimum_defect_reduction':min(r['defect_reduction'] for r in checks if r['mode']=='on')}
    result['passed']=result['relative_error']<=.001 and result['base_current_relative_change']<=1e-8 and result['adjoint_relative_residual']<=1e-10
    p.audit.write(OUT/'fixed_state_result.json',result);print(json.dumps(result,indent=2),flush=True)


def dc():
    p.audit.verify(FREEZE)
    assert p.audit.read(OUT/'fixed_state_result.json')['passed']
    rows=[]
    for name,amp in p.CASES:
        dest=LOCAL/'dc'/name;source=p.LOCAL/name/'charge.csv'
        s=execute(dest/'config.json',True,source)
        cfg=p.audit.read(dest/'config.json');cfg.pop('output_state_file')
        cfg['state_file']=str(dest/'state.csv');cfg['solver']['carrier_row_convergence']['mode']='report'
        cfg['solver']['global_continuity_closure']=p.audit.read(p.CONTRACT)['global_audit_profile']
        cfg.update(simulation_type='newton_carrier_term_probe',output_csv=str(dest/'terms.csv'),carrier_term_probe={'solved_equation_terms':True})
        a=p.execute(cfg,dest/'acceptance.json',source)
        cc=s['contact_currents_A_per_um'];id=cc['drain'];kcl=abs(math.fsum(cc.values()))/abs(id)
        old=next(r for r in p.audit.rows(p.OUT/'case_ledger.csv') if r['name']==name)
        r={'name':name,'amplitude_cm_3':amp,'current_A_per_um':id,'exit_code':s['exit_code'],'iterations':s['iterations'],'reason':s['convergence_reason'],
            'local_violations':a['carrier_row_convergence']['violation_count'],'kcl_over_Id':kcl,'current_relative_change':id/float(old['current_A_per_um'])-1,
            'electron_global_qualified':a['global_continuity_closure']['electron']['qualified'],'hole_global_qualified':a['global_continuity_closure']['hole']['qualified'],
            'qualified':s['converged'] and a['carrier_row_convergence']['satisfied'] and a['global_continuity_closure']['satisfied'] and kcl<=1e-8}
        rows.append(r);print(r,flush=True)
    p.audit.write_csv(OUT/'dc_ledger.csv',rows)
    by={r['name']:r['current_A_per_um'] for r in rows};pred=p.audit.read(OUT/'fixed_state_result.json')['prediction_A_per_um_at_1e13']
    cal=[]
    for suffix,amp in (('1e13',1e13),('5e12',5e12)):
        response=(by['plus_'+suffix]-by['minus_'+suffix])/2
        cal.append({'amplitude_cm_3':amp,'central_response_A_per_um':response,'prediction_A_per_um':pred*amp/1e13,'relative_error':abs(response/(pred*amp/1e13)-1)})
    p.audit.write_csv(OUT/'dc_calibration.csv',cal)
    p.audit.write(OUT/'result.json',{'strict_states':sum(r['qualified'] for r in rows),'max_fd_adjoint_relative':max(r['relative_error'] for r in cal),
        'max_current_relative_change':max(abs(r['current_relative_change']) for r in rows),
        'passed':all(r['qualified'] and abs(r['current_relative_change'])<=1e-8 for r in rows) and all(r['relative_error']<=.001 for r in cal),
        'new_nonlinear_solves':5,'new_sentaurus_runs':0,'m82_released':False,'m83_released':False})


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=('build','prepare','check','dc','verify'))
    {'build':build,'prepare':prepare,'check':check,'dc':dc,'verify':lambda:p.audit.verify(FREEZE)}[parser.parse_args().action]()
