"""Frozen conservative dielectric candidate and interface-only causal control.

No production source/default is edited. Native geometry retains Vela epsilon0;
charge constants, charge/SRH volumes and native effective mobility are fixed.
"""
import argparse, copy, math, subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
import localize_simplemos_remaining_poisson_20260907 as l
import check_simplemos_joint_geometry_20260907 as jc

a=l.a;d=l.d;v=l.v;q=l.q;p=l.p
LOCAL=l.LOCAL/'dielectric';OUT=l.OUT/'dielectric';RUNNER=LOCAL/'runner.exe'
AXES=('all_native_K','Si_SiO2_only')

HEADER=r'''
#include <cstdlib>
#include <fstream>
#include <vector>
#include <cmath>
#include <stdexcept>
namespace vela_dielectric_candidate {
inline double factor(std::size_t edge, std::size_t size) {
    const char* value=std::getenv("VELA_CANDIDATE_DIELECTRIC_ALPHA");
    if (!value || std::stod(value)==0.) return 1.;
    const double alpha=std::stod(value);
    if (!std::isfinite(alpha)) throw std::runtime_error("Nonfinite dielectric amplitude");
    static const auto ratios=[&] {
        const char* path=std::getenv("VELA_CANDIDATE_DIELECTRIC_RATIOS");
        if (!path) throw std::runtime_error("Missing dielectric ratios");
        std::ifstream in(path);std::vector<double> rows(size);
        for(std::size_t k=0;k<size;++k) {
            std::size_t id;
            if (!(in>>id>>rows[k]) || id!=k || !std::isfinite(rows[k]) || rows[k]<0.)
                throw std::runtime_error("Invalid dielectric ratio row");
        }
        std::string extra;if(in>>extra)throw std::runtime_error("Extra dielectric rows");
        return rows;
    }();
    if(size!=ratios.size())throw std::runtime_error("Dielectric mesh mismatch");
    const double result=1.+alpha*(ratios.at(edge)-1.);
    if(!(result>=0.) || !std::isfinite(result))throw std::runtime_error("Invalid dielectric multiplier");
    return result;
}
}
'''

def build():
    a.verify(l.OUT/'ledger_evidence.json');LOCAL.mkdir(parents=True,exist_ok=False)
    source=(q.LOCAL/'CoupledDDAssembler.cpp').read_text()
    marker='        const auto appendStencilNode = [&](Index node) {';assert source.count(marker)==1
    source=HEADER+source.replace(marker,'        kernel.poissonCoupling *= vela_dielectric_candidate::factor(edgeId, mesh_.numEdges());\n\n'+marker)
    (LOCAL/'CoupledDDAssembler.cpp').write_text(source,newline='\n')
    cmd=copy.deepcopy(a.read(q.LOCAL/'compile_command.json'));cmd[cmd.index('-o')+1]=str(LOCAL/'CoupledDDAssembler.o');cmd[cmd.index('-c')+1]=str(LOCAL/'CoupledDDAssembler.cpp')
    a.write(LOCAL/'compile_command.json',cmd)
    r=subprocess.run(cmd,cwd=p.REPO/'build-release',env=v.environment(),capture_output=True,text=True)
    (LOCAL/'build.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr[-3000:]
    cmd=a.read(q.LOCAL/'link_command.json');cmd=[str(LOCAL/'CoupledDDAssembler.o') if x==str(q.LOCAL/'CoupledDDAssembler.o') else x for x in cmd];cmd[-1]=str(RUNNER)
    a.write(LOCAL/'link_command.json',cmd);r=subprocess.run(cmd,env=v.environment(),capture_output=True,text=True)
    (LOCAL/'link.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr[-3000:]
    print('Isolated dielectric runner built',flush=True)

def env(c,axis,alpha):
    e=l.env(c,1);e.update(VELA_CANDIDATE_DIELECTRIC_ALPHA=format(alpha,'.17g'),
                         VELA_CANDIDATE_DIELECTRIC_RATIOS=str(LOCAL/'ratios'/(c['device']+'_'+axis+'.txt')))
    return e

def execute(path,c,axis,alpha):return v.execute(path,RUNNER,env(c,axis,alpha))

def prepare():
    a.verify(l.OUT/'ledger_evidence.json');assert all(x['qualified']=='True' for x in a.rows(l.OUT/'ledger_checks.csv'))
    cases=a.read(l.OUT/'ledger_contract.json')['cases'];files=[Path(__file__).resolve(),l.OUT/'ledger_evidence.json',RUNNER,q.RUNNER,p.REPO/'build-release/libvela_core.a'];geotests=[]
    for device in ('n19','n23'):
        geo,xy,el,vol,K,parts,*_=l.g.geometry(device)
        c=next(c for c in cases if c['device']==device);edges=a.rows(l.LOCAL/c['key']/'native_mu/edges.csv')
        for axis in AXES:
            values=[];source=np.zeros(geo.count);absolute=0.;changed=0
            for e in edges:
                i,j=sorted((int(e['node0']),int(e['node1'])));old=float(-geo.matrices['legacy'][i,j]);target=float(-K[i,j])
                selected=axis=='all_native_K' or {x['material'] for x in parts[(i,j)]}=={'Si','SiO2'}
                new=target if selected else old
                if old==0:assert abs(new)<=1e-25;ratio=1.
                else:ratio=new/old
                assert math.isfinite(ratio) and ratio>=0
                assert abs(old*ratio-new)<=1e-14*max(abs(new),abs(old),1e-300)
                values.append(f"{e['edge_id']} {ratio:.17g}\n")
                change=new-old;changed+=abs(change)>1e-25
                value=change*math.sin(.41*(i+1)) # independent arbitrary scalar edge flux
                source[i]+=value;source[j]-=value;absolute+=2*abs(value)
            path=LOCAL/'ratios'/(device+'_'+axis+'.txt');path.parent.mkdir(parents=True,exist_ok=True);path.write_text(''.join(values),newline='\n');files.append(path)
            conservation=abs(math.fsum(source))/max(absolute,1e-300);assert conservation<=1e-14
            geotests.append(dict(device=device,axis=axis,edges=len(edges),changed_edges=changed,edge_pair_conservation_relative=conservation,all_multipliers_nonnegative=True,qualified=True))
    for c in cases:
        cfg=v.config(c);root=LOCAL/c['key'];state=q.LOCAL/c['key']/'replacement/state.csv';c['dielectric_probes']=[];c['dielectric_jobs']=[]
        for label,axis,alpha in [('zero',AXES[0],0.)]+[(axis,axis,1.) for axis in AXES]:
            dest=root/'fixed'/label;paths=v.probes(cfg,dest,state)
            deck=copy.deepcopy(a.read(dest/'edges.json'));deck.update(simulation_type='joint_geometry_jvp',output_csv=str(dest/'jvp.csv'));a.write(dest/'jvp.json',deck);paths.append(dest/'jvp.json')
            c['dielectric_probes'] += [dict(path=str(x),axis=axis,alpha=alpha) for x in paths if x.stem in ('functional','edges','jvp') or (x.stem=='adjoint' and label=='zero')];files+=paths
        jobs=[('zero',AXES[0],0.,state)]
        for axis in AXES:
            jobs += [(axis+'/'+label,axis,alpha,state) for label,alpha in [('plus_full',.001),('minus_full',-.001),('plus_half',.0005),('minus_half',-.0005),('replacement',1.)]]
            jobs += [(axis+'/replacement_native',axis,1.,v.m.previous.LOCAL/c['key']/'native_referenced_joint/state.csv')]
        for label,axis,alpha,initial in jobs:
            dest=root/label;deck=copy.deepcopy(cfg);deck.update(state_file=str(initial),output_state_file=str(dest/'state.csv'))
            a.write(dest/'config.json',deck);v.post_config(deck,dest)
            c['dielectric_jobs'].append(dict(label=label,axis=axis,alpha=alpha,dest=str(dest)))
            files += [dest/'config.json',dest/'all_row.json',dest/'acceptance_edges.json',initial]
            if 'replacement' in label:files+=v.probes(cfg,dest/'post',dest/'state.csv')
    a.write_csv(OUT/'geometry_properties.csv',geotests)
    a.write(OUT/'contract.json',dict(cases=cases,gates=a.read(v.OUT/'contract.json')['gates'],small_DC=72,finite_DC=32,
        candidate='K(alpha)=K_legacy+alpha*(K_native_geometry-K_legacy), using unchanged Vela epsilon0. Shared cached symmetric edge coefficient in residual and analytic Jacobian.',
        control='Si_SiO2_only selects the 20 shared Si/SiO2 interface edges; all_native_K changes all native dielectric edge geometry. Each edge is conservative.',
        fixed='Qualified joint geometry and native effective electron/hole mobility, all three signed Si Poisson charge volumes, original SRH volume, Vela q/kb/epsilon0, boundaries and acceptance.',
        promotion='Numerical qualification and dual initialization required; compare full result to both mobility-only and production joint baselines, high/low NWell pairing. Interface-only arm is a localization control, not a production model.',
        source_gate='Poisson differential relative <=1e-8 on all noncontact nodes; carrier rows and current unchanged bitwise at fixed state.',
        defaults_changed=False,acceptance_changed=False))
    files += [OUT/'contract.json',OUT/'geometry_properties.csv']+[x for x in LOCAL.glob('*') if x.is_file()]
    d.matrix.freeze(OUT/'freeze.json',files);print('Frozen: 72 small DC and 32 gated finite DC, two independent directions',flush=True)

def preflight():
    a.verify(OUT/'freeze.json');cases=a.read(OUT/'contract.json')['cases']
    def one(c):
        for j in c['dielectric_probes']:
            s=execute(Path(j['path']),c,j['axis'],j['alpha']);assert s['exit_code']==0,(j,s)
        print('Dielectric probes',c['key'],flush=True)
    with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(one,cases))
    factor=a.read(d.matrix.p.CONTRACT)['frozen_poisson_residual_per_physical_charge'];checks=[];preds=[];jvp=[]
    for c in cases:
        geo,xy,el,vol,K,parts,*_=l.g.geometry(c['device']);root=LOCAL/c['key']/'fixed'
        zeroR=d.ordered(root/'zero/residual.csv',geo.count);zeroE=a.rows(root/'zero/edges.csv');state=l.live(zeroE,geo.count)
        adj=d.ordered(root/'zero/adjoint.csv',geo.count);lam=np.array([float(x['lambda_poisson']) for x in adj]);lam[geo.contact_nodes]=0
        zjp,jbase=jc.jvp_metrics(root/'zero/jvp.csv',c['key'],'zero');jvp+=zjp
        old=a.read(q.LOCAL/c['key']/'replacement/config.status.json')['contact_currents_A_per_um']['drain']
        s0=a.read(root/'zero/functional.status.json');drift=abs(s0['current_A_per_um']/old-1)
        for axis in AXES:
            ratios={int(x.split()[0]):float(x.split()[1]) for x in (LOCAL/'ratios'/(c['device']+'_'+axis+'.txt')).read_text().splitlines()}
            source=np.zeros(geo.count)
            for e in zeroE:
                i,j=int(e['node0']),int(e['node1']);oldK=float(-geo.matrices['legacy'][i,j]);newK=float(-K[i,j])
                selected=axis=='all_native_K' or {x['material'] for x in parts[tuple(sorted((i,j)))]}=={'Si','SiO2'}
                delta=(newK-oldK) if selected else 0.;value=delta*(state[0,i]-state[0,j])*factor
                source[i]+=value;source[j]-=value
                assert abs(oldK*(ratios[int(e['edge_id'])]-1)-delta)<=1e-14*max(abs(oldK),1e-300)
            source[geo.contact_nodes]=0
            u=d.ordered(root/axis/'residual.csv',geo.count);ue=a.rows(root/axis/'edges.csv')
            actual=np.array([float(x['psi_residual'])-float(y['psi_residual']) for x,y in zip(u,zeroR)])
            err=float(np.linalg.norm((actual-source)[geo.free])/np.linalg.norm(source[geo.free]))
            unchanged=all(x[k]==y[k] for x,y in zip(u,zeroR) for k in ('phin_residual','phip_residual'))
            edge_unchanged=all(x[k]==y[k] for x,y in zip(ue,zeroE) for k in ('electron_flux','hole_flux','couple_m','electron_mobility_m2_V_s','hole_mobility_m2_V_s'))
            s=a.read(root/axis/'functional.status.json');port=abs(s['current_A_per_um']/s['contact_current_extractor_A_per_um']-1)
            unchanged_current=s['current_A_per_um']==s0['current_A_per_um']
            jp,_=jc.jvp_metrics(root/axis/'jvp.csv',c['key'],axis,jbase);jvp+=jp
            prediction=-math.fsum(float(x*y) for x,y in zip(lam,source))
            preds.append(dict(key=c['key'],axis=axis,unit_prediction_A_per_um=prediction,unit_direct_A_per_um=0.,relative_to_Id=prediction/old))
            checks.append(dict(key=c['key'],axis=axis,Poisson_source_relative=err,carrier_rows_bit_identical=unchanged,transport_bit_identical=edge_unchanged,
                fixed_current_bit_identical=unchanged_current,port_relative=port,zero_operator_Id_relative=drift,
                qualified=err<=1e-8 and unchanged and edge_unchanged and unchanged_current and port<=1e-8 and drift<=1e-8 and all(x['qualified'] for x in zjp+jp)))
    for name,rows in [('preflight',checks),('predictions',preds),('jvp',jvp)]:a.write_csv(OUT/(name+'.csv'),rows)
    d.matrix.freeze(OUT/'preflight_evidence.json',[OUT/'freeze.json',OUT/'preflight.csv',OUT/'predictions.csv',OUT/'jvp.csv']+[x for c in cases for x in (LOCAL/c['key']/'fixed').rglob('*') if x.is_file()])
    print('Dielectric preflight qualified',sum(x['qualified'] for x in checks),'/',len(checks),flush=True)

def qualify(c,dest):
    """Pure independent acceptance; runner/env chosen explicitly before this call."""
    s=a.read(dest/'config.status.json');ps=a.read(dest/'all_row.status.json');es=a.read(dest/'acceptance_edges.status.json')
    geo,mask=v.m.previous.prior.support(c);terms=d.ordered(dest/'all_row.csv',geo.count);edges=a.rows(dest/'acceptance_edges.csv');contacts=set(map(int,geo.contact_nodes))
    ratios=[];closure={}
    for r,keep in zip(terms,mask):
        if keep:
            for car in ('electron','hole'):
                scale=max(float(r[car+'_flux_abs_sum']),abs(float(r[car+'_recombination'])),abs(float(r[car+'_impact'])))
                ratios.append(abs(float(r[car+'_residual']))/scale if scale else math.inf)
    assert len(ratios)==1814
    for car in ('electron','hole'):
        flux=math.fsum((int(int(x['node0']) in contacts)-int(int(x['node1']) in contacts))*float(x[car+'_flux']) for x in edges)
        source=math.fsum(float(r[car+'_recombination'])+float(r[car+'_impact']) for i,r in enumerate(terms) if i not in contacts)
        ratio=abs(flux-source)/max(abs(flux),abs(source),1e-10)
        closure[car]=dict(contact_flux=flux,integrated_source=source,ratio=ratio,qualified=abs(source)>=1e-10,satisfied=abs(source)<1e-10 or ratio<=1e-6)
    cc=s.get('contact_currents_A_per_um',{});Id=cc.get('drain',math.nan);kcl=abs(math.fsum(cc.values()))/abs(Id);bad=sum(x>1e-6 for x in ratios)
    result=dict(qualified=s['exit_code']==ps['exit_code']==es['exit_code']==0 and s.get('converged',False) and bad==0 and all(x['satisfied'] for x in closure.values()) and kcl<=1e-8,
        current_A_per_um=Id,max_row_ratio=max(ratios),row_violations=bad,kcl_over_Id=kcl,iterations=s.get('iterations'),failure=s.get('failure_reason',''),
        global_electron_ratio=closure['electron']['ratio'],global_hole_ratio=closure['hole']['ratio'])
    record=dict(result=result,closure=closure);target=dest/'independent_acceptance.json'
    if target.exists():assert a.read(target)==record
    else:a.write(target,record)
    return result

def run(finite=False):
    a.verify(OUT/'freeze.json');a.verify(OUT/'preflight_evidence.json');assert all(x['qualified']=='True' for x in a.rows(OUT/'preflight.csv'))
    if finite:
        a.verify(OUT/'response_evidence.json');assert all(x['qualified']=='True' for x in a.rows(OUT/'response.csv'))
    def one(c):
        rows=[]
        for j in c['dielectric_jobs']:
            if ('replacement' in j['label'])!=finite:continue
            dest=Path(j['dest']);s=execute(dest/'config.json',c,j['axis'],j['alpha'])
            row=dict(key=c['key'],axis=j['axis'],label=j['label'],qualified=False,failure=s.get('failure_reason',''))
            if (dest/'state.csv').exists():
                for name in ('all_row','acceptance_edges'):execute(dest/(name+'.json'),c,j['axis'],j['alpha'])
                row.update(qualify(c,dest))
                if finite:
                    for name in ('functional','edges','terms'):execute(dest/'post'/(name+'.json'),c,j['axis'],j['alpha'])
            rows.append(row);print(c['key'],j['label'],row['qualified'],row.get('max_row_ratio'),flush=True)
        return rows
    with ThreadPoolExecutor(max_workers=2) as pool:rows=[r for group in pool.map(one,a.read(OUT/'contract.json')['cases']) for r in group]
    keys=list(dict.fromkeys(k for r in rows for k in r));a.write_csv(OUT/('finite_dc.csv' if finite else 'small_dc.csv'),[{k:r.get(k,'') for k in keys} for r in rows])

def response():
    dc=a.rows(OUT/'small_dc.csv');rows=[]
    for c in a.read(OUT/'contract.json')['cases']:
        data={r['label']:r for r in dc if r['key']==c['key']};center=float(data['zero']['current_A_per_um']);old=a.read(q.LOCAL/c['key']/'replacement/config.status.json')['contact_currents_A_per_um']['drain'];drift=abs(center-old)
        for axis in AXES:
            pred=float(next(r['unit_prediction_A_per_um'] for r in a.rows(OUT/'predictions.csv') if r['key']==c['key'] and r['axis']==axis))
            val={k:float(data[axis+'/'+k]['current_A_per_um']) for k in ('plus_full','minus_full','plus_half','minus_half')}
            odd={k:(val['plus_'+k]-val['minus_'+k])/2 for k in ('full','half')};lin=abs(odd['full']/(2*odd['half'])-1)
            for amp,alpha in [('full',.001),('half',.0005)]:
                error=abs(odd[amp]/(alpha*pred)-1);even=abs((val['plus_'+amp]+val['minus_'+amp])/2-center)/abs(odd[amp]);snr=abs(odd[amp])/max(drift,1e-300);dd=abs(math.log10(center/old))
                signs=(val['plus_'+amp]-center)*pred>0 and (val['minus_'+amp]-center)*pred<0
                passed=data['zero']['qualified']=='True' and all(data[axis+'/'+k+'_'+amp]['qualified']=='True' for k in ('plus','minus')) and error<=.001 and lin<=.001 and even<=.01 and snr>=100 and dd<=1e-5 and signs
                rows.append(dict(key=c['key'],axis=axis,amplitude=amp,prediction_relative_error=error,two_amplitude_relative=lin,even_over_odd=even,signal_over_zero_drift=snr,zero_drift_dex=dd,unit_prediction_A_per_um=pred,actual_odd_A_per_um=odd[amp],qualified=passed))
    a.write_csv(OUT/'response.csv',rows)
    d.matrix.freeze(OUT/'response_evidence.json',[OUT/'freeze.json',OUT/'preflight_evidence.json',OUT/'small_dc.csv',OUT/'response.csv']+[x for c in a.read(OUT/'contract.json')['cases'] for j in c['dielectric_jobs'] if 'replacement' not in j['label'] for x in Path(j['dest']).rglob('*') if x.is_file()])
    print('Dielectric same-source response',sum(x['qualified'] for x in rows),'/',len(rows),flush=True)

def analyze():
    a.verify(OUT/'response_evidence.json');dc=a.rows(OUT/'finite_dc.csv');comparison=[];dual=[];local=[];pairing=[];ports=[]
    for c in a.read(OUT/'contract.json')['cases']:
        geo,mask=v.m.previous.prior.support(c);meta={k:c[k] for k in ('key','device','vg','vd')}
        oldstate=d.ordered(q.LOCAL/c['key']/'replacement/state.csv',geo.count)
        oldterms=d.ordered(q.LOCAL/c['key']/'replacement/all_row.csv',geo.count)
        oldfields={int(r['node_id']):r for r in a.rows(q.OUT/'local_fields.csv') if r['key']==c['key']}
        muId=a.read(q.LOCAL/c['key']/'replacement/config.status.json')['contact_currents_A_per_um']['drain'];prodId=a.read(v.LOCAL/c['key']/'from_vela/config.status.json')['contact_currents_A_per_um']['drain']
        for axis in AXES:
            root=LOCAL/c['key']/axis;states=[d.ordered(root/name/'state.csv',geo.count) for name in ('replacement','replacement_native')];infos=[qualify(c,root/name) for name in ('replacement','replacement_native')]
            diff=v.m.previous.prior.delta_states(*states,mask);Id=infos[0]['current_A_per_um'];err=abs(Id/infos[1]['current_A_per_um']-1)
            ok=all(x['qualified'] for x in infos) and err<=1e-6 and max(diff[k] for k in ('psi_max_V','phin_max_V','phip_max_V'))<=1e-6 and diff['density_max_relative']<=1e-4
            dual.append(dict(**meta,axis=axis,Id_relative=err,**diff,qualified=ok))
            pred=float(next(r['unit_prediction_A_per_um'] for r in a.rows(OUT/'predictions.csv') if r['key']==c['key'] and r['axis']==axis))
            comparison.append(dict(**meta,axis=axis,native_Id_A_per_um=c['native_Id_A_per_um'],production_Id_A_per_um=prodId,mobility_Id_A_per_um=muId,candidate_Id_A_per_um=Id,
                production_relative_error=prodId/c['native_Id_A_per_um']-1,mobility_relative_error=muId/c['native_Id_A_per_um']-1,candidate_relative_error=Id/c['native_Id_A_per_um']-1,
                predicted_finite_delta_A_per_um=pred,actual_finite_delta_A_per_um=Id-muId,finite_to_tangent_relative=(Id-muId)/pred-1,qualified=ok))
            newterms=d.ordered(root/'replacement/all_row.csv',geo.count)
            for node in (1000,1009):
                r=dict(**meta,axis=axis,node_id=node,qualified=ok);old=oldfields[node]
                for f in ('psi','phin','phip'):
                    change=float(v.m.previous.prior.physical(states[0][node],f)-v.m.previous.prior.physical(oldstate[node],f))
                    r[f+'_old_delta_V']=float(old[f+'_candidate_delta_V']);r[f+'_new_delta_V']=r[f+'_old_delta_V']+change
                for car,f in [('electron','electrons_m3'),('hole','holes_m3')]:
                    native_density=float(oldstate[node][f])/(1+float(old[car+'_candidate_density_relative']))
                    r[car+'_new_density_relative']=float(states[0][node][f])/native_density-1
                rate=float(old['candidate_SRH_cm3_s'])*float(newterms[node]['electron_recombination'])/float(oldterms[node]['electron_recombination'])
                r.update(new_SRH_cm3_s=rate,native_plot_SRH_cm3_s=float(old['native_plot_SRH_cm3_s']),old_SRH_relative=float(old['candidate_SRH_relative']),new_SRH_relative=rate/float(old['native_plot_SRH_cm3_s'])-1);local.append(r)
            for name in ('replacement','replacement_native'):
                status=a.read(root/name/'post/functional.status.json');error=abs(status['current_A_per_um']/status['contact_current_extractor_A_per_um']-1)
                ports.append(dict(**meta,axis=axis,initialization=name,relative_error=error,qualified=status['exit_code']==0 and error<=1e-8))
    for axis in AXES:
        for vg in (.8,1.):
            for vd in (.05,1.):
                pair={dev:next(r for r in comparison if r['device']==dev and r['vg']==vg and r['vd']==vd and r['axis']==axis) for dev in ('n19','n23')}
                ratio=pair['n23']['native_Id_A_per_um']/pair['n19']['native_Id_A_per_um']
                new=math.log10((pair['n23']['candidate_Id_A_per_um']/pair['n19']['candidate_Id_A_per_um'])/ratio)
                for base in ('production','mobility'):
                    old=math.log10((pair['n23'][base+'_Id_A_per_um']/pair['n19'][base+'_Id_A_per_um'])/ratio)
                    hi=abs(pair['n23']['candidate_relative_error'])<abs(pair['n23'][base+'_relative_error']);lo=abs(pair['n19']['candidate_relative_error'])<=abs(pair['n19'][base+'_relative_error'])
                    pairing.append(dict(axis=axis,vg=vg,vd=vd,baseline=base,old_pair_error_dex=old,new_pair_error_dex=new,high_NWell_improved=hi,low_NWell_not_worse=lo,pair_improved=abs(new)<abs(old),qualified=hi and lo and abs(new)<abs(old)))
    for name,rows in [('comparison',comparison),('dual_init',dual),('local_fields',local),('pairing',pairing),('ports',ports)]:a.write_csv(OUT/(name+'.csv'),rows)
    a.write(OUT/'summary.json',dict(finite_qualified=sum(r['qualified']=='True' for r in dc),finite_count=len(dc),dual_qualified=sum(x['qualified'] for x in dual),port_qualified=sum(x['qualified'] for x in ports),
        all_native_K_promoted=all(x['qualified'] for x in pairing if x['axis']=='all_native_K') and all(x['qualified'] for x in dual+ports if x['axis']=='all_native_K'),
        production_changed=False,full_curve_started=False))
    d.matrix.freeze(OUT/'validation_evidence.json',[OUT/'freeze.json',OUT/'response_evidence.json',OUT/'summary.json']+list(OUT.glob('*.csv'))+[x for c in a.read(OUT/'contract.json')['cases'] for j in c['dielectric_jobs'] if 'replacement' in j['label'] for x in Path(j['dest']).rglob('*') if x.is_file()])
    print(a.read(OUT/'summary.json'),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('build','prepare','preflight','small','response','finite','analyze'));action=parser.parse_args().action
    if action in ('small','finite'):run(action=='finite')
    else:globals()[action]()
