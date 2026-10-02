"""Recalibrate the same conservative channel-flux direction on qualified Masetti states."""
import argparse
import copy
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import json
import math
import subprocess
import prepare_simplemos_masetti_local_20260907 as p
import validate_simplemos_local_conservative_flux as old

a=p.a;d=p.d;LOCAL=p.LOCAL/'channel_response';OUT=p.OUT/'channel_response';RUNNER=LOCAL/'runner.exe'


def env(root,edge=None,flux=0.):
    e=p.prev.env('masetti_refined',root)
    e.pop(old.EDGE_ENV,None);e.pop(old.FLUX_ENV,None)
    if edge is not None:e[old.EDGE_ENV]=str(edge);e[old.FLUX_ENV]=format(flux,'.17g')
    return e


def build():
    a.verify(p.prev.OUT/'validation_evidence.json');LOCAL.mkdir(parents=True,exist_ok=False)
    args=a.read(p.prev.p.prior.LOCAL/'build_command.json')
    src=Path(next(x for x in args if Path(x).name=='CoupledDDAssembler.cpp'))
    text=src.read_text()
    for marker in ('            r(phinOffset() + i) += nFlux;', '            terms[static_cast<Index>(i)].electronFlux += nFlux;'):
        assert text.count(marker)==1;text=text.replace(marker,'            nFlux += diagnosticConstantElectronEdgeFlux(e);\n'+marker)
    (LOCAL/'CoupledDDAssembler.cpp').write_text(old.HELPER+text,newline='\n')
    # Reuse the already independently checked contact implementation of exactly
    # this prescribed conservative source. Its production flux remains unchanged.
    contact=old.LOCAL/'ContactCurrent.cpp';assert contact.exists()
    args=[str(LOCAL/'CoupledDDAssembler.cpp') if x==str(src) else x for x in args]
    args.insert(next(i for i,x in enumerate(args) if x.endswith('libvela_core.a')),str(contact));args[-1]=str(RUNNER)
    a.write(LOCAL/'build_command.json',args)
    done=subprocess.run(args,env=env(LOCAL),capture_output=True,text=True)
    (LOCAL/'build.log').write_text(done.stdout+done.stderr);assert done.returncode==0,done.stderr[-3000:]
    print('Built isolated conservative-source runner with existing four-refinement Newton',flush=True)


def prepare():
    a.verify(p.OUT/'vela_freeze.json');files=[Path(__file__).resolve(),RUNNER,p.OUT/'vela_freeze.json',old.LOCAL/'ContactCurrent.cpp'];cases=[]
    for c in a.read(p.OUT/'vela_contract.json')['cases']:
        if c['vg']!=1.:continue
        root=LOCAL/c['key'];base=Path(c['base']);deck=a.read(base/'config.json');deck['state_file']=str(base/'state.csv')
        es=a.rows(p.LOCAL/'vela'/c['key']/'strict/edges.csv');edge=next(x for x in es if (int(x['node0']),int(x['node1']))==(320,324))
        adj=d.ordered(p.LOCAL/'vela'/c['key']/'adjoint.csv',len(a.rows(base/'state.csv')))
        assert not adj[320]['contacts'] and not adj[324]['contacts']
        current=a.read(base/'result.json')['Id_A_per_um'];amount=.0005*current
        conversion=d.fixed.Q*1e-6*float(edge['electron_particle_line_flux_per_m_s'])/float(edge['electron_flux'])
        scaled=amount/conversion;native=old.native_flux_for_current(amount)
        prediction=-(float(adj[320]['lambda_electron'])-float(adj[324]['lambda_electron']))*scaled
        jobs=[]
        for label,mult in [('zero',0.),('plus_full',1.),('minus_full',-1.),('plus_half',.5),('minus_half',-.5)]:
            dest=root/label;cfg=copy.deepcopy(deck);cfg['output_state_file']=str(dest/'state.csv')
            cfg['solver']['local_update_diagnostics']['csv_file']=str(dest/'updates.csv');a.write(dest/'config.json',cfg)
            q=copy.deepcopy(cfg);q.pop('output_state_file');q['solver'].pop('local_update_diagnostics')
            q.update(simulation_type='newton_carrier_term_probe',state_file=str(dest/'state.csv'),output_csv=str(dest/'all_row.csv'),carrier_term_probe={'solved_equation_terms':True})
            q['solver']['carrier_row_convergence']['mode']='report';q['solver']['global_continuity_closure']=dict(mode='enforce',tolerance=1e-6,source_floor=1e-10)
            a.write(dest/'all_row.json',q);jobs.append(dict(label=label,mult=mult,config=str(dest/'config.json')))
            files += [dest/'config.json',dest/'all_row.json']
        for label,mult in [('zero',0.),('full',1.)]:
            cfg=copy.deepcopy(deck);cfg.pop('output_state_file');cfg['solver'].pop('local_update_diagnostics')
            cfg.update(simulation_type='terminal_current_functional_probe',contact='drain',residual_output_csv=str(root/'preflight'/label/'residual.csv'))
            a.write(root/'preflight'/label/'config.json',cfg);files.append(root/'preflight'/label/'config.json')
        cases.append(dict(key=c['key'],device=c['device'],vd=c['vd'],vg=c['vg'],base=str(base),edge=int(edge['edge_id']),node0=320,node1=324,
            base_Id_A_per_um=current,native_amplitude=native,scaled_amplitude=scaled,current_amplitude_A_per_um=amount,prediction_A_per_um=prediction,jobs=jobs))
        files += [base/'config.json',base/'state.csv',base/'result.json',p.LOCAL/'vela'/c['key']/'adjoint.csv',p.LOCAL/'vela'/c['key']/'strict/edges.csv']
    assert len(cases)==4
    a.write(OUT/'contract.json',dict(status='frozen_before_execution',cases=cases,new_self_consistent_runs=20,
        scope='n19/n23 x Vd=.05/1 at Vg=1. Same free-node channel edge 320-324. Recalibrate fresh Masetti weights, not reuse PhuMob weights.',
        source='State-independent conservative electron edge flux, + at node0 and - at node1 before contact rows. Same source in post row diagnostics; no direct contact term on this free edge.',
        amplitude='Full/half = .05%/.025% of each frozen Vela Id; fixed before execution.',
        gates=dict(source_relative=1e-8,fd_relative=.001,two_amplitude_relative=.001,even_fraction=.01,minimum_signal_to_zero_drift=100,zero_drift_dex=1e-5,kcl_over_Id=1e-8),
        original_all_row_and_global_gates_unchanged=True,production_changes=False))
    files.append(OUT/'contract.json');files += [Path(x) for x in a.read(LOCAL/'build_command.json') if Path(x).is_file() and Path(x).is_relative_to(p.REPO)]
    d.matrix.freeze(OUT/'freeze.json',files);print('Frozen 20 Masetti channel response solves',flush=True)


def execute(path,edge,flux):
    target=path.with_suffix('.status.json')
    if target.exists():return a.read(target)
    done=subprocess.run([str(RUNNER),'--config',str(path),'--log','off'],env=env(path.parent,edge,flux),capture_output=True,text=True)
    path.with_suffix('.stdout.txt').write_text(done.stdout);path.with_suffix('.stderr.txt').write_text(done.stderr)
    s=json.loads(done.stdout.strip().splitlines()[-1]);s['exit_code']=done.returncode;a.write(target,s);return s


def run():
    a.verify(OUT/'freeze.json')
    def one(c):
        root=LOCAL/c['key'];probes={}
        for label,mult in [('zero',0.),('full',1.)]:
            s=execute(root/'preflight'/label/'config.json',c['edge'],mult*c['native_amplitude']);assert s['exit_code']==0
            probes[label]=s
        count=len(a.rows(Path(c['base'])/'state.csv'));z=d.ordered(root/'preflight/zero/residual.csv',count);f=d.ordered(root/'preflight/full/residual.csv',count)
        oldz=d.ordered(p.LOCAL/'vela'/c['key']/'strict/residual.csv',count)
        assert all(x[k]==y[k] for x,y in zip(z,oldz) for k in ('psi_residual','phin_residual','phip_residual'))
        for i,(x,y) in enumerate(zip(z,f)):
            assert x['psi_residual']==y['psi_residual'] and x['phip_residual']==y['phip_residual']
            expect=c['scaled_amplitude']*(int(i==320)-int(i==324));actual=float(y['phin_residual'])-float(x['phin_residual'])
            assert abs(actual-expect)<=1e-8*abs(c['scaled_amplitude'])
        assert probes['zero']['current_A_per_um']==probes['full']['current_A_per_um']
        results=[]
        for job in c['jobs']:
            path=Path(job['config']);flux=job['mult']*c['native_amplitude'];s=execute(path,c['edge'],flux);q=execute(path.parent/'all_row.json',c['edge'],flux)
            cc=s['contact_currents_A_per_um'];kcl=abs(math.fsum(cc.values()))/abs(cc['drain']);g=q['carrier_row_convergence']
            geo=d.matrix.spatial.m73.Geometry(c['device'])
            mask=d.matrix.spatial.old.m78.supports(c['device'],geo,.05)[0]['all_si'].copy();mask[geo.contact_nodes]=False
            terms=d.ordered(path.parent/'all_row.csv',geo.count)
            positive_scales=all(max(float(r[k+'_flux_abs_sum']),abs(float(r[k+'_recombination'])),abs(float(r[k+'_impact'])))>0 for r,keep in zip(terms,mask) if keep for k in ('electron','hole'))
            assert 2*int(sum(mask))==1814
            qualified=s['exit_code']==q['exit_code']==0 and s['converged'] and g['satisfied'] and g['qualified_row_count']==1814 and positive_scales and q['global_continuity_closure']['satisfied'] and kcl<=1e-8
            results.append(dict(key=c['key'],label=job['label'],current_A_per_um=cc['drain'],converged=s['converged'],iterations=s['iterations'],failure=s['failure_reason'],qualified=qualified,row_violations=g['violation_count'],max_row_ratio=g['max_ratio'],kcl_over_Id=kcl))
        print(c['key'],'qualified',sum(x['qualified'] for x in results),'/ 5',flush=True);return results
    with ThreadPoolExecutor(max_workers=2) as pool:rows=[x for group in pool.map(one,a.read(OUT/'contract.json')['cases']) for x in group]
    a.write_csv(OUT/'dc.csv',rows)


def analyze():
    a.verify(OUT/'freeze.json');rows=a.rows(OUT/'dc.csv');result=[]
    for c in a.read(OUT/'contract.json')['cases']:
        values={x['label']:x for x in rows if x['key']==c['key']};zero=float(values['zero']['current_A_per_um'])
        full=(float(values['plus_full']['current_A_per_um'])-float(values['minus_full']['current_A_per_um']))/2
        half=(float(values['plus_half']['current_A_per_um'])-float(values['minus_half']['current_A_per_um']))/2
        linearity=abs(full/(2*half)-1)
        drift=abs(zero-c['base_Id_A_per_um']);driftdex=abs(math.log10(zero/c['base_Id_A_per_um']))
        for name,scale in [('full',1.),('half',.5)]:
            plus=float(values['plus_'+name]['current_A_per_um']);minus=float(values['minus_'+name]['current_A_per_um']);odd=(plus-minus)/2;even=(plus+minus)/2-zero
            error=abs(odd/(scale*c['prediction_A_per_um'])-1);ratio=abs(odd)/max(drift,1e-300)
            sign=(plus-zero)*c['prediction_A_per_um']>0 and (minus-zero)*c['prediction_A_per_um']<0
            passed=all(x['qualified']=='True' for x in values.values()) and error<=.001 and linearity<=.001 and abs(even/odd)<=.01 and ratio>=100 and driftdex<=1e-5 and sign
            result.append(dict(key=c['key'],amplitude=name,prediction_A_per_um=scale*c['prediction_A_per_um'],odd_A_per_um=odd,
                prediction_relative_error=error,two_amplitude_relative=linearity,even_over_odd=abs(even/odd),signal_over_zero_drift=ratio,zero_drift_dex=driftdex,signs_correct=sign,qualified=passed))
    a.write_csv(OUT/'calibration.csv',result);print('Qualified channel response amplitudes',sum(x['qualified'] for x in result),'/',len(result),flush=True)


if __name__=='__main__':
    p0=argparse.ArgumentParser();p0.add_argument('action',choices=('build','prepare','run','analyze'));globals()[p0.parse_args().action]()
