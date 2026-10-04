"""Isolated constant electron-edge flux: direct terminal term plus DD feedback."""
import argparse
import copy
import json
import math
from pathlib import Path
import subprocess
import numpy as np
import decompose_simplemos_calibrated_transport as prior

a=prior.a
m=prior.matrix
REPO=prior.REPO
LOCAL=REPO/'build-release/simplemos_local_conservative_flux_20260906'
OUT=prior.ROOT/'local_conservative_flux'
RUNNER=LOCAL/'local_flux_runner.exe'
CONTRACT=OUT/'contract.json'
FREEZE=OUT/'freeze.json'
DOC=REPO/'docs/validation/simplemos_local_conservative_flux_2026-09-06.md'
EDGE_ENV='VELA_DIAGNOSTIC_CONSTANT_ELECTRON_EDGE'
FLUX_ENV='VELA_DIAGNOSTIC_CONSTANT_ELECTRON_FLUX_NATIVE'
# PhysicalUnitSystem::tcadInternal: current-density scale 1e4, length scale 1e-6.
LINE_FACTOR=1e4*1e-6
HELPER=r'''
#include <cstdlib>
#include <cmath>
#include <stdexcept>
namespace {
double diagnosticConstantElectronEdgeFlux(std::size_t edge) {
    const char* e=std::getenv("VELA_DIAGNOSTIC_CONSTANT_ELECTRON_EDGE");
    const char* f=std::getenv("VELA_DIAGNOSTIC_CONSTANT_ELECTRON_FLUX_NATIVE");
    if (!e && !f) return 0.0;
    if (!e || !f) throw std::runtime_error("Incomplete diagnostic edge flux");
    const auto id=std::stoll(e);
    const double value=std::stod(f);
    if (id<0 || !std::isfinite(value)) throw std::runtime_error("Invalid diagnostic edge flux");
    return static_cast<std::size_t>(id)==edge ? value : 0.0;
}
}
'''


def env(edge=None,flux=0.):
    result=m.chain.env(True)
    result['VELA_LINEAR_SOLVER']='sparselu'
    result.pop(EDGE_ENV,None);result.pop(FLUX_ENV,None)
    if edge is not None:
        result[EDGE_ENV]=str(edge);result[FLUX_ENV]=format(flux,'.17g')
    return result


def build():
    a.verify(prior.FREEZE);a.verify(prior.OUT/'evidence.json')
    if RUNNER.exists():raise FileExistsError(RUNNER)
    LOCAL.mkdir(parents=True,exist_ok=True)
    source=(m.chain.LOCAL/'CoupledDDAssembler.cpp').read_text()
    for marker in ('            r(phinOffset() + i) += nFlux;', '            terms[static_cast<Index>(i)].electronFlux += nFlux;'):
        assert source.count(marker)==1
        source=source.replace(marker,'            nFlux += diagnosticConstantElectronEdgeFlux(e);\n'+marker)
    (LOCAL/'CoupledDDAssembler.cpp').write_text(HELPER+source,newline='\n')
    contact=(REPO/'src/post/ContactCurrent.cpp').read_text()
    marker='        // Algebraic SG split: J = J_drift + J_diffusion.'
    assert contact.count(marker)==1
    contact=contact.replace(marker,r'''
        // Same prescribed conservative line flux as the residual. A contact
        // receives the explicit direct term, so residual feedback alone is not Id.
        const double diagnosticDensityFlux = diagnosticConstantElectronEdgeFlux(e) / couple_[e];
        electronContinuityFlux01 += diagnosticDensityFlux;
        electronContinuityFluxLongDouble01 += diagnosticDensityFlux;
        electronFlux01 -= diagnosticDensityFlux;
'''+marker)
    (LOCAL/'ContactCurrent.cpp').write_text(HELPER+contact,newline='\n')
    # Reuse the frozen read-only acceptance extension; it does not alter Newton.
    auditdir=REPO/'build-release/simplemos_convergence_audit'
    runner=(m.chain.LOCAL/'runner.cpp').read_text()
    old=(auditdir/'runner.cpp').read_text()
    start='nlohmann::json runNewtonCarrierTermProbe('
    end='void writeSgEdgeFluxProbeCsv'
    runner=runner[:runner.index(start)]+old[old.index(start):old.index(end)]+runner[runner.index(end):]
    (LOCAL/'runner.cpp').write_text(runner,newline='\n')
    args=a.read(m.chain.LOCAL/'build_command.json')
    args=[str(LOCAL/'runner.cpp') if Path(v).name=='runner.cpp' else str(LOCAL/'CoupledDDAssembler.cpp') if Path(v).name=='CoupledDDAssembler.cpp' else v for v in args]
    args.insert(1,'-I'+str(auditdir/'include'))
    index=next(i for i,v in enumerate(args) if v.endswith('libvela_core.a'))
    args[index:index]=[str(LOCAL/'ContactCurrent.cpp'),str(auditdir/'NewtonSolver.cpp')]
    args[-1]=str(RUNNER)
    a.write(LOCAL/'build_command.json',args)
    proc=subprocess.run(args,env=env(),capture_output=True,text=True)
    (LOCAL/'build.log').write_text(proc.stdout+proc.stderr)
    if proc.returncode:raise RuntimeError(proc.stderr[-3000:])
    print('Built isolated conservative flux runner; production files unchanged',flush=True)


def prediction(lam_i,lam_j,scaled_flux,direct):
    feedback=-(lam_i-lam_j)*scaled_flux
    return {'feedback_A_per_um':feedback,'direct_A_per_um':direct,'total_A_per_um':feedback+direct}


def native_flux_for_current(current_A_per_um):
    return current_A_per_um/(prior.fixed.Q*1e-6*LINE_FACTOR)


def prepare():
    a.verify(prior.FREEZE);a.verify(prior.OUT/'evidence.json')
    cases=[];files=[Path(__file__).resolve(),RUNNER,prior.OUT/'evidence.json',
        REPO/'tests/regression/test_simplemos_local_conservative_flux.py',REPO/'src/post/ContactCurrent.cpp']
    for item in a.read(LOCAL/'build_command.json'):
        if Path(item).is_file():files.append(Path(item))
    files.append(REPO/'build-release/simplemos_convergence_audit/include/vela/solver/NewtonSolver.h')
    for c0 in a.read(prior.CONTRACT)['cases']:
        if c0['device']!='n23':continue
        c=copy.deepcopy(c0);key=c['case'];root=LOCAL/key
        base=m.spatial.precision.LOCAL/key
        deck=a.read(base/'config.json');deck['state_file']=str(base/'state.csv')
        adj=prior.ordered(Path(c['adjoint']),c['nodes'])
        edges={int(r['edge_id']):r for r in a.rows(prior.LOCAL/key/'strict/edges.csv')}
        regions=[]
        for name,edge,relative in (('drain_contact',2392,.05),('channel',965,.0005)):
            e=edges[edge];i,j=int(e['node0']),int(e['node1'])
            # Physical current per normalized residual unit, independently available
            # from production edge flux and particle line flux at the same state.
            conv=prior.fixed.Q*1e-6*float(e['electron_particle_line_flux_per_m_s'])/float(e['electron_flux'])
            amount=relative*float(c['current_A_per_um'])
            residual_amount=amount/conv
            li=0. if adj[i]['contacts'] else float(adj[i]['lambda_electron'])
            lj=0. if adj[j]['contacts'] else float(adj[j]['lambda_electron'])
            direct=-amount*(int('drain' in adj[i]['contacts'].split(';'))-int('drain' in adj[j]['contacts'].split(';')))
            pred=prediction(li,lj,residual_amount,direct)
            # Recover native unscaled line flux from the diagnostic physical ratio.
            # This is C0*D0 times the normalized residual amplitude.
            scale=conv/(prior.fixed.Q*1e-6*LINE_FACTOR)
            native=native_flux_for_current(amount)
            assert abs(native/residual_amount/scale-1)<1e-12
            r=dict(name=name,edge=edge,node0=i,node1=j,relative_current=relative,current_amplitude_A_per_um=amount,
                scaled_amplitude=residual_amount,native_amplitude=native,conversion_A_per_um=conv,**pred)
            regions.append(r)
            for label,mult in (('plus_full',1.),('minus_full',-1.),('plus_half',.5),('minus_half',-.5)):
                dest=root/name/label
                d=copy.deepcopy(deck);d['output_state_file']=str(dest/'state.csv')
                a.write(dest/'config.json',d);files.append(dest/'config.json')
        c['regions']=regions
        d=copy.deepcopy(deck);d['output_state_file']=str(root/'zero/state.csv')
        a.write(root/'zero/config.json',d);files.append(root/'zero/config.json')
        files += [base/'config.json',base/'state.csv',Path(c['adjoint']),prior.LOCAL/key/'strict/edges.csv']
        files += [Path(deck[k]) for k in ('mesh_file','materials_file','node_doping_file')]
        cases.append(c)
    a.write(CONTRACT,{'status':'frozen_before_execution','cases':cases,
        'source':'Constant electron particle line flux on one existing Si edge, + at node0 and - at node1 before Dirichlet replacement. Same flux in carrier diagnostics and conventional contact current; Poisson/hole/Jacobian unchanged.',
        'amplitudes':'Before execution: contact edge 5%/2.5% of strict Id; channel edge .05%/.025%. Independent frozen amplitudes accommodate direct/feedback cancellation, not fitted FD outcomes.',
        'native_line_flux_to_particles_per_m_s':LINE_FACTOR,
        'gates':{'fd_relative':.001,'two_amplitude_relative':.001,'even_fraction':.01,'minimum_signal_to_zero_drift':100,'zero_drift_dex':1e-5,'kcl_over_Id':1e-8,'source_relative':1e-8},
        'global_profile':a.read(m.CONTRACT)['global_audit_profile'],
        'new_nonlinear_solves':18,'new_sentaurus_runs':0,'production_changes':False,'m82_released':False,'m83_released':False})
    files.append(CONTRACT)
    a.write(FREEZE,{'input_hashes':{(a.rel(f) if f.is_relative_to(REPO) else str(f)):a.sha(f) for f in sorted(set(files))}})
    print('Frozen 18 DC states and two conservative edge directions at both Vd',flush=True)


def freeze_prepared():
    """Finish an interrupted pre-execution manifest without rewriting decks."""
    assert CONTRACT.exists() and not FREEZE.exists()
    cfg=a.read(CONTRACT)
    decks=list(LOCAL.rglob('config.json'))
    assert len(decks)==18 and not list(LOCAL.rglob('*.status.json'))
    files=[Path(__file__).resolve(),RUNNER,CONTRACT,prior.OUT/'evidence.json',
        REPO/'tests/regression/test_simplemos_local_conservative_flux.py',REPO/'src/post/ContactCurrent.cpp',
        REPO/'build-release/simplemos_convergence_audit/include/vela/solver/NewtonSolver.h']+decks
    files += [Path(v) for v in a.read(LOCAL/'build_command.json') if Path(v).is_file()]
    for path in decks:
        deck=a.read(path)
        files += [Path(deck[k]) for k in ('state_file','mesh_file','materials_file','node_doping_file')]
    for c in cfg['cases']:
        files += [Path(c['adjoint']),m.spatial.precision.LOCAL/c['case']/'config.json',prior.LOCAL/c['case']/'strict/edges.csv']
    a.write(FREEZE,{'input_hashes':{(a.rel(f) if f.is_relative_to(REPO) else str(f)):a.sha(f) for f in sorted(set(files))}})
    print('Completed pre-execution manifest with absolute paths for external toolchain inputs',flush=True)


def execute(path,edge=None,native=0.):
    statuspath=path.with_suffix('.status.json')
    if statuspath.exists():return a.read(statuspath)
    proc=subprocess.run([str(RUNNER),'--config',str(path),'--log','off'],env=env(edge,native),capture_output=True,text=True)
    path.with_suffix('.stdout.txt').write_text(proc.stdout)
    path.with_suffix('.stderr.txt').write_text(proc.stderr)
    s=json.loads(proc.stdout.strip().splitlines()[-1]);s['exit_code']=proc.returncode
    a.write(statuspath,s);return s


def preflight():
    a.verify(FREEZE);checks=[]
    for c in a.read(CONTRACT)['cases']:
        root=LOCAL/c['case']
        deck=a.read(root/'zero/config.json');deck.pop('output_state_file')
        outputs={}
        for name,edge,native in [('zero',None,0.)]+[(r['name'],r['edge'],r['native_amplitude']) for r in c['regions']]:
            dest=root/'preflight'/name
            d=dict(deck,simulation_type='newton_residual_probe',output_csv=str(dest/'residual.csv'))
            a.write(dest/'residual.json',d);s=execute(dest/'residual.json',edge,native);assert s['exit_code']==0
            outputs[name]=prior.ordered(dest/'residual.csv',c['nodes'])
            d=dict(deck,simulation_type='terminal_current_functional_probe',contact='drain')
            a.write(dest/'functional.json',d);s=execute(dest/'functional.json',edge,native);assert s['exit_code']==0
        # Fresh zero replay against the sealed unmodified residual.
        z=outputs['zero'];ref=prior.ordered(prior.LOCAL/c['case']/'strict/residual.csv',c['nodes'])
        assert all(u[k]==v[k] for u,v in zip(z,ref) for k in ('psi_residual','phin_residual','phip_residual'))
        for r in c['regions']:
            rows=outputs[r['name']]
            actual=np.array([float(u['phin_residual'])-float(v['phin_residual']) for u,v in zip(rows,z)])
            expected=np.zeros(c['nodes']);expected[r['node0']]=r['scaled_amplitude'];expected[r['node1']]=-r['scaled_amplitude']
            adj=prior.ordered(Path(c['adjoint']),c['nodes'])
            expected[[i for i,v in enumerate(adj) if v['contacts']]]=0
            error=float(np.linalg.norm(actual-expected)/np.linalg.norm(expected))
            assert error<=1e-8
            assert all(u[k]==v[k] for u,v in zip(rows,z) for k in ('psi_residual','phip_residual','net_doping_m3','ni_eff_m3'))
            before=a.read(root/'preflight/zero/functional.status.json')['current_A_per_um']
            after=a.read(root/'preflight'/r['name']/'functional.status.json')['current_A_per_um']
            assert abs(after-before-r['direct_A_per_um'])<=1e-8*r['current_amplitude_A_per_um']
            checks.append({'case':c['case'],'region':r['name'],'source_relative':error,'direct_current_A_per_um':after-before,'passed':True})
        print(c['case'],'residual source and direct current preflight passed',flush=True)
    a.write_csv(OUT/'preflight.csv',checks)


def run():
    a.verify(FREEZE);assert len(a.rows(OUT/'preflight.csv'))==4
    ledger=[]
    for c in a.read(CONTRACT)['cases']:
        root=LOCAL/c['case']
        jobs=[('zero','zero',None,0.)]
        for r in c['regions']:
            jobs += [(r['name'],label,r['edge'],r['native_amplitude']*mult) for label,mult in (('plus_full',1.),('minus_full',-1.),('plus_half',.5),('minus_half',-.5))]
        for region,label,edge,native in jobs:
            dest=root/'zero' if region=='zero' else root/region/label
            s=execute(dest/'config.json',edge,native)
            d=a.read(dest/'config.json');d.pop('output_state_file');d['state_file']=str(dest/'state.csv')
            d['solver']['carrier_row_convergence']['mode']='report'
            d['solver']['global_continuity_closure']=a.read(CONTRACT)['global_profile']
            d.update(simulation_type='newton_carrier_term_probe',output_csv=str(dest/'terms.csv'),carrier_term_probe={'solved_equation_terms':True})
            if not (dest/'acceptance.json').exists():a.write(dest/'acceptance.json',d)
            audit=execute(dest/'acceptance.json',edge,native)
            cc=s['contact_currents_A_per_um'];current=cc['drain'];kcl=abs(math.fsum(cc.values()))/abs(current)
            local,glob=audit['carrier_row_convergence'],audit['global_continuity_closure']
            r={'case':c['case'],'region':region,'label':label,'current_A_per_um':current,'iterations':s['iterations'],'reason':s['convergence_reason'],
                'exit_code':s['exit_code'],'local_violations':local['violation_count'],'kcl_over_Id':kcl,
                'electron_global_qualified':glob['electron']['qualified'],'hole_global_qualified':glob['hole']['qualified'],
                'qualified':s['exit_code']==0 and audit['exit_code']==0 and s['converged'] and local['satisfied'] and glob['satisfied'] and kcl<=1e-8}
            ledger.append(r);print(json.dumps(r),flush=True)
    a.write_csv(OUT/'dc.csv',ledger)


def analyze():
    a.verify(FREEZE);ledger=a.rows(OUT/'dc.csv');cal=[]
    for c in a.read(CONTRACT)['cases']:
        group=[r for r in ledger if r['case']==c['case']]
        z=next(r for r in group if r['region']=='zero');zero=float(z['current_A_per_um'])
        drift=abs(zero-float(c['current_A_per_um']))
        for region in c['regions']:
            rows=[r for r in group if r['region']==region['name']]
            vals={r['label']:float(r['current_A_per_um']) for r in rows}
            f=(vals['plus_full']-vals['minus_full'])/2;h=(vals['plus_half']-vals['minus_half'])/2
            linearity=abs(2*h/f-1) if f else float('inf')
            qualified=len(rows)==4 and z['qualified']=='True' and all(r['qualified']=='True' for r in rows)
            for suffix,mult,fd in (('full',1.,f),('half',.5,h)):
                target=region['total_A_per_um']*mult
                err=abs(fd/target-1)
                even=abs((vals['plus_'+suffix]+vals['minus_'+suffix])/2-zero)/abs(fd)
                signal=abs(fd)/max(drift,1e-300)
                signs=(vals['plus_'+suffix]-zero)*target>0 and (vals['minus_'+suffix]-zero)*target<0
                passed=qualified and err<=.001 and linearity<=.001 and even<=.01 and signal>=100 and signs and abs(math.log10(zero/float(c['current_A_per_um'])))<=1e-5
                cal.append({'case':c['case'],'region':region['name'],'amplitude':suffix,'predicted_A_per_um':target,'actual_A_per_um':fd,
                    'relative_error':err,'two_amplitude_relative':linearity,'even_fraction':even,'signal_to_zero_drift':signal,
                    'response_over_injected_current':fd/(region['current_amplitude_A_per_um']*mult),'passed':passed})
    a.write_csv(OUT/'calibration.csv',cal)
    result={'strict_states':sum(r['qualified']=='True' for r in ledger),'passed_amplitudes':sum(r['passed'] for r in cal),'amplitudes':len(cal),
        'new_nonlinear_solves':len(ledger),'new_sentaurus_runs':0,'m82_released':False,'m83_released':False}
    a.write(OUT/'result.json',result);print(json.dumps(result,indent=2),flush=True)


def seal():
    a.verify(FREEZE)
    m.freeze(OUT/'evidence.json',[DOC]+[f for root in (LOCAL,OUT) for f in root.rglob('*') if f.is_file()])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=('build','prepare','freeze_prepared','preflight','run','analyze','seal','verify'))
    action=parser.parse_args().action
    if action=='verify':a.verify(FREEZE);a.verify(OUT/'evidence.json');print('Local flux evidence verified')
    else:globals()[action]()
