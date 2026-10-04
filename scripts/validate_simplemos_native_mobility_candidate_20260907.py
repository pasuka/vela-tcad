"""Same-direction adjoint/FD calibration followed by gated finite replacement."""
import argparse,copy,math
from concurrent.futures import ThreadPoolExecutor
from decimal import Decimal
from pathlib import Path
import numpy as np
import build_simplemos_native_mobility_candidate_20260907 as b
import check_simplemos_joint_geometry_20260907 as jc
import audit_simplemos_masetti_box_mobility_20260907 as native
import complete_simplemos_production_migration_20260907 as v
a=v.a;d=v.d;p=v.p;LOCAL=b.LOCAL;OUT=b.OUT

def env(c,alpha):
    e=v.environment();e.update(VELA_CANDIDATE_MOBILITY_ALPHA=format(alpha,'.17g'),VELA_CANDIDATE_MOBILITY_RATIOS=str(LOCAL/'ratios'/(c['device']+'.txt')));return e

def execute(path,c,alpha):return v.execute(path,b.RUNNER,env(c,alpha))

def prepare():
    a.verify(v.OUT/'validation_evidence.json');assert a.read(v.OUT/'summary.json')['production_qualified']
    cases=a.read(v.OUT/'contract.json')['cases'];geometry=[];files=[Path(__file__).resolve(),Path(b.__file__).resolve(),v.OUT/'validation_evidence.json',b.RUNNER]
    for device in ('n19','n23'):
        geo,weights,checks=native.geometry(device);assert all(x['qualified'] for x in checks)
        _,_,_,_,_,parts,_,_=native.geom.geometry(device)
        c=next(c for c in cases if c['device']==device)
        edges=a.rows(v.LOCAL/c['key']/'fixed/joint/edges.csv');lines=[]
        for e in edges:
            i,j=int(e['node0']),int(e['node1']);g=math.fsum(x['coefficient'] for x in parts[(i,j)] if x['material']=='Si');ratios=[]
            for car,k in [('electron','e'),('hole','h')]:
                mu=float(e[car+'_mobility_m2_V_s'])*1e4
                ratio=weights[(i,j)][k]/g/mu if g>0 and mu>0 else 1.
                assert ratio>0 and math.isfinite(ratio);ratios.append(ratio)
                geometry.append(dict(device=device,edge_id=int(e['edge_id']),carrier=car,baseline_cm2_Vs=mu,native_effective_cm2_Vs=mu*ratio,ratio=ratio,native_Si_g=g))
            lines.append(f"{e['edge_id']} {ratios[0]:.17g} {ratios[1]:.17g}\n")
        path=LOCAL/'ratios'/(device+'.txt');path.parent.mkdir(parents=True,exist_ok=True);path.write_text(''.join(lines),newline='\n');files.append(path)
    for c in cases:
        c['probes']=[];c['jobs']=[];cfg=v.config(c);root=LOCAL/c['key'];state=v.LOCAL/c['key']/'from_vela/state.csv'
        # Calibrate weights at the actual production-qualified joint state.
        for label,alpha in [('zero',0.),('unit',1.)]:
            dest=root/'fixed'/label;paths=v.probes(cfg,dest,state)
            j=copy.deepcopy(a.read(dest/'edges.json'));j.update(simulation_type='joint_geometry_jvp',output_csv=str(dest/'jvp.csv'));a.write(dest/'jvp.json',j);paths.append(dest/'jvp.json')
            c['probes'] += [dict(path=str(x),alpha=alpha) for x in paths];files+=paths
        for label,alpha,initial in [('zero',0.,state),('plus_full',.001,state),('minus_full',-.001,state),('plus_half',.0005,state),('minus_half',-.0005,state),
                                    ('replacement',1.,state),('replacement_native',1.,v.m.previous.LOCAL/c['key']/'native_referenced_joint/state.csv')]:
            dest=root/label;deck=copy.deepcopy(cfg);deck.update(state_file=str(initial),output_state_file=str(dest/'state.csv'))
            a.write(dest/'config.json',deck);v.post_config(deck,dest)
            c['jobs'].append(dict(label=label,alpha=alpha,dest=str(dest)));files += [dest/'config.json',dest/'all_row.json',initial]
            if label.startswith('replacement'):
                files+=v.probes(cfg,dest/'post',dest/'state.csv')
    a.write_csv(OUT/'mobility_ratios.csv',geometry)
    a.write(OUT/'contract.json',dict(status='frozen_before_execution',cases=cases,small_DC=40,finite_DC=16,gates=a.read(v.OUT/'contract.json')['gates'],
        candidate='Fixed electron and hole native Si box-effective mobility; interpolation mu(alpha)=mu_Vela*(1+alpha*(mu_native/mu_Vela-1)). Only this constant coefficient changes.',
        geometry='Qualified joint transport and Poisson charge geometry. Original dielectric K and SRH source volumes.',
        promotion='First pass all eight same-source +/- two-amplitude response gates. Then finite replacement with dual initialization, high-NWell absolute improvement, low-NWell not worse, NWell pairing improvement at both Vg and Vd.',
        acceptance_changes=False,production_mobility_default_changed=False))
    files += [OUT/'contract.json',OUT/'mobility_ratios.csv']+[x for x in LOCAL.rglob('*') if x.is_file() and x.suffix in ('.cpp','.h','.json','.log')]
    d.matrix.freeze(OUT/'freeze.json',files);print('Frozen native mobility candidate: 40 small-perturbation and 16 gated replacement DC',flush=True)

def preflight():
    a.verify(OUT/'freeze.json');cases=a.read(OUT/'contract.json')['cases']
    def one(c):
        for j in c['probes']:
            s=execute(Path(j['path']),c,j['alpha']);assert s['exit_code']==0,(j,s)
        print('Candidate fixed-state probes',c['key'],flush=True)
    with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(one,cases))
    checks=[];predictions=[];jvp=[]
    for c in cases:
        root=LOCAL/c['key']/'fixed';geo,mask=v.m.previous.prior.support(c)
        zeroE=a.rows(root/'zero/edges.csv');unitE=a.rows(root/'unit/edges.csv');zeroR=d.ordered(root/'zero/residual.csv',geo.count);unitR=d.ordered(root/'unit/residual.csv',geo.count)
        ratios={int(x.split()[0]):[float(y)-1 for y in x.split()[1:]] for x in (LOCAL/'ratios'/(c['device']+'.txt')).read_text().splitlines()}
        adj=d.ordered(root/'zero/adjoint.csv',geo.count);lam=np.array([[float(x[k]) for x in adj] for k in ('lambda_poisson','lambda_electron','lambda_hole')]);lam[:,geo.contact_nodes]=0
        source=np.zeros((3,geo.count));flow=np.zeros_like(source);absolute=np.zeros_like(source);direct=[]
        cfg=a.read(root/'zero/edges.json');mesh=a.read(Path(cfg['mesh_file']));drain=set(next(x['node_ids'] for x in mesh['contacts'] if x['name']=='drain'))
        edge_error=0.
        for e,u in zip(zeroE,unitE):
            i,j=int(e['node0']),int(e['node1']);rr=ratios[int(e['edge_id'])]
            assert e['couple_m']==u['couple_m']
            for block,car in [(1,'electron'),(2,'hole')]:
                delta=rr[block-1];f=float(e[car+'_flux']);s=f*delta
                source[block,i]+=s;source[block,j]-=s;flow[block,i]+=abs(s);flow[block,j]+=abs(s);absolute[block,i]+=abs(f);absolute[block,j]+=abs(f)
                expected=f*(1+delta);edge_error=max(edge_error,abs(float(u[car+'_flux'])-expected)/max(abs(expected),1e-300))
            direct.append(-d.fixed.Q*1e-6*(int(i in drain)-int(j in drain))*(float(e['electron_particle_line_flux_per_m_s'])*rr[0]-float(e['hole_particle_line_flux_per_m_s'])*rr[1]))
        maxbound=0.
        for i,(z,u) in enumerate(zip(zeroR,unitR)):
            assert z['psi_residual']==u['psi_residual']
            if not mask[i]:continue
            for block,column in [(1,'phin_residual'),(2,'phip_residual')]:
                actual=float(u[column])-float(z[column]);bound=max(1e-8*flow[block,i]+64*np.finfo(float).eps*absolute[block,i],1e-25)
                maxbound=max(maxbound,abs(actual-source[block,i])/bound)
        source[:,geo.contact_nodes]=0;feedback=-math.fsum(float(x*y) for x,y in zip(lam.flat,source.flat));di=math.fsum(direct)
        predictions.append(dict(key=c['key'],unit_prediction_A_per_um=feedback+di,unit_feedback_A_per_um=feedback,unit_direct_A_per_um=di))
        ports=[]
        for label in ('zero','unit'):
            s=a.read(root/label/'functional.status.json');ports.append(abs(s['current_A_per_um']/s['contact_current_extractor_A_per_um']-1))
        zjp,jbase=jc.jvp_metrics(root/'zero/jvp.csv',c['key'],'zero');ujp,_=jc.jvp_metrics(root/'unit/jvp.csv',c['key'],'unit',jbase);jvp+=zjp+ujp
        # Source-free zero must match the actual production operator, irrespective of compile unit layout.
        prodId=a.read(v.LOCAL/c['key']/'from_vela/config.status.json')['contact_currents_A_per_um']['drain']
        zeroId=a.read(root/'zero/functional.status.json')['current_A_per_um'];drift=abs(zeroId/prodId-1)
        checks.append(dict(key=c['key'],source_error_over_bound=maxbound,edge_relative_error=edge_error,port_relative=max(ports),zero_operator_Id_relative=drift,
            qualified=maxbound<=1 and edge_error<=1e-10 and max(ports)<=1e-8 and drift<=1e-8 and all(x['qualified'] for x in zjp+ujp)))
    a.write_csv(OUT/'preflight.csv',checks);a.write_csv(OUT/'predictions.csv',predictions);a.write_csv(OUT/'jvp.csv',jvp)
    d.matrix.freeze(OUT/'preflight_evidence.json',[OUT/'freeze.json',OUT/'preflight.csv',OUT/'predictions.csv',OUT/'jvp.csv']+[x for c in cases for x in (LOCAL/c['key']/'fixed').rglob('*') if x.is_file()])
    print('Candidate preflight passed:',all(x['qualified'] for x in checks),flush=True)

def run(finite=False):
    a.verify(OUT/'freeze.json');a.verify(OUT/'preflight_evidence.json');assert all(x['qualified']=='True' for x in a.rows(OUT/'preflight.csv'))
    if finite:
        a.verify(OUT/'response_evidence.json');assert all(x['qualified']=='True' for x in a.rows(OUT/'response.csv'))
    def one(c):
        rows=[]
        for j in c['jobs']:
            if j['label'].startswith('replacement')!=finite:continue
            dest=Path(j['dest']);s=execute(dest/'config.json',c,j['alpha']);row=dict(key=c['key'],device=c['device'],vg=c['vg'],vd=c['vd'],label=j['label'],qualified=False,failure=s.get('failure_reason',''))
            if (dest/'state.csv').exists():
                execute(dest/'all_row.json',c,j['alpha']);row.update(v.qualify(c,dest))
                if finite:
                    for name in ('functional','edges','terms','mobility'):execute(dest/'post'/(name+'.json'),c,j['alpha'])
            rows.append(row);print(c['key'],j['label'],row['qualified'],row.get('max_row_ratio'),flush=True)
        return rows
    with ThreadPoolExecutor(max_workers=2) as pool:rows=[r for group in pool.map(one,a.read(OUT/'contract.json')['cases']) for r in group]
    keys=list(dict.fromkeys(k for r in rows for k in r));a.write_csv(OUT/('finite_dc.csv' if finite else 'small_dc.csv'),[{k:r.get(k,'') for k in keys} for r in rows])

def response():
    dc=a.rows(OUT/'small_dc.csv');rows=[]
    for c in a.read(OUT/'contract.json')['cases']:
        data={r['label']:r for r in dc if r['key']==c['key']};values={k:float(r['current_A_per_um']) for k,r in data.items()}
        pred=float(next(r['unit_prediction_A_per_um'] for r in a.rows(OUT/'predictions.csv') if r['key']==c['key']))
        center=values['zero'];old=a.read(v.LOCAL/c['key']/'from_vela/config.status.json')['contact_currents_A_per_um']['drain'];drift=abs(center-old)
        odd={k:(values['plus_'+k]-values['minus_'+k])/2 for k in ('full','half')};lin=abs(odd['full']/(2*odd['half'])-1)
        for k,alpha in [('full',.001),('half',.0005)]:
            error=abs(odd[k]/(alpha*pred)-1);even=abs((values['plus_'+k]+values['minus_'+k])/2-center)/abs(odd[k]);snr=abs(odd[k])/max(drift,1e-300);dd=abs(math.log10(center/old))
            signs=(values['plus_'+k]-center)*pred>0 and (values['minus_'+k]-center)*pred<0
            ok=all(data[t]['qualified']=='True' for t in ('zero','plus_'+k,'minus_'+k)) and error<=.001 and lin<=.001 and even<=.01 and snr>=100 and dd<=1e-5 and signs
            rows.append(dict(key=c['key'],amplitude=k,prediction_relative_error=error,two_amplitude_relative=lin,even_over_odd=even,signal_over_zero_drift=snr,zero_drift_dex=dd,signs_correct=signs,unit_prediction_A_per_um=pred,actual_odd_A_per_um=odd[k],qualified=ok))
    a.write_csv(OUT/'response.csv',rows);d.matrix.freeze(OUT/'response_evidence.json',[OUT/'freeze.json',OUT/'preflight_evidence.json',OUT/'response.csv',OUT/'small_dc.csv']+[x for c in a.read(OUT/'contract.json')['cases'] for j in c['jobs'] if not j['label'].startswith('replacement') for x in Path(j['dest']).rglob('*') if x.is_file()])
    print('Mobility same-source amplitudes passed:',sum(x['qualified'] for x in rows),'/',len(rows),flush=True)

def analyze():
    a.verify(OUT/'response_evidence.json');dc=a.rows(OUT/'finite_dc.csv');comparison=[];dual=[];fields=[];pairing=[]
    for c in a.read(OUT/'contract.json')['cases']:
        geo,mask=v.m.previous.prior.support(c);root=LOCAL/c['key'];meta={k:c[k] for k in ('key','device','vg','vd')}
        initial=d.ordered(v.LOCAL/c['key']/'from_vela/state.csv',geo.count)
        s1=d.ordered(root/'replacement/state.csv',geo.count);s2=d.ordered(root/'replacement_native/state.csv',geo.count)
        diff=v.m.previous.prior.delta_states(s1,s2,mask);infos=[v.qualify(c,root/name) for name in ('replacement','replacement_native')];Id=infos[0]['current_A_per_um']
        err=abs(Id/infos[1]['current_A_per_um']-1);ok=all(x['qualified'] for x in infos) and err<=1e-6 and max(diff[k] for k in ('psi_max_V','phin_max_V','phip_max_V'))<=1e-6 and diff['density_max_relative']<=1e-4
        dual.append(dict(**meta,Id_relative=err,**diff,qualified=ok))
        base=a.read(v.LOCAL/c['key']/'from_vela/config.status.json')['contact_currents_A_per_um']['drain']
        comparison.append(dict(**meta,native_Id_A_per_um=c['native_Id_A_per_um'],base_Id_A_per_um=base,candidate_Id_A_per_um=Id,base_relative_error=base/c['native_Id_A_per_um']-1,candidate_relative_error=Id/c['native_Id_A_per_um']-1,qualified=ok))
        previous={int(r['node_id']):r for r in a.rows(v.m.previous.OUT/'local_node_ledger.csv') if r['key']==c['key']}
        baseTerms=d.ordered(v.LOCAL/c['key']/'from_vela/all_row.csv',geo.count);newTerms=d.ordered(root/'replacement/all_row.csv',geo.count)
        for i in (1000,1009):
            old=previous[i];row=dict(**meta,node_id=i,qualified_state=infos[0]['qualified'])
            for f in ('psi','phin','phip'):
                change=float(v.m.previous.prior.physical(s1[i],f)-v.m.previous.prior.physical(initial[i],f))
                # Native value recovered from the earlier immutable local ledger/state.
                oldState=d.ordered(Path(c['root'])/'joint/replacement/state.csv',geo.count)[i]
                nativeValue=v.m.previous.prior.physical(oldState,f)-Decimal(old[f+'_delta_V'])
                row[f+'_candidate_delta_V']=float(v.m.previous.prior.physical(s1[i],f)-nativeValue);row[f+'_change_V']=change
            for car,f in [('electron','electrons_m3'),('hole','holes_m3')]:row[car+'_candidate_density_relative']=float(s1[i][f])/float(old[car+'_native_density_m3'])-1
            # Same SRH volume/scaling: scale the old physical rate by the solved equation term ratio.
            oldTerm=d.ordered(v.m.previous.LOCAL/c['key']/'joint/terms.csv',geo.count)[i]
            rate=float(old['Vela_SRH_cm3_s'])*float(newTerms[i]['electron_recombination'])/float(oldTerm['electron_recombination'])
            row.update(candidate_SRH_cm3_s=rate,native_plot_SRH_cm3_s=float(old['native_plot_SRH_cm3_s']),candidate_SRH_relative=rate/float(old['native_plot_SRH_cm3_s'])-1)
            fields.append(row)
    for vg in (.8,1.):
        for vd in (.05,1.):
            pair={dev:next(r for r in comparison if r['device']==dev and r['vg']==vg and r['vd']==vd) for dev in ('n19','n23')}
            nativeRatio=pair['n23']['native_Id_A_per_um']/pair['n19']['native_Id_A_per_um']
            old=math.log10((pair['n23']['base_Id_A_per_um']/pair['n19']['base_Id_A_per_um'])/nativeRatio);new=math.log10((pair['n23']['candidate_Id_A_per_um']/pair['n19']['candidate_Id_A_per_um'])/nativeRatio)
            hi=abs(pair['n23']['candidate_relative_error'])<abs(pair['n23']['base_relative_error']);lo=abs(pair['n19']['candidate_relative_error'])<=abs(pair['n19']['base_relative_error'])
            pairing.append(dict(vg=vg,vd=vd,base_pair_error_dex=old,candidate_pair_error_dex=new,high_NWell_improved=hi,low_NWell_not_worse=lo,pair_improved=abs(new)<abs(old),qualified=hi and lo and abs(new)<abs(old)))
    for name,data in [('comparison',comparison),('dual_init',dual),('local_fields',fields),('NWell_pairing',pairing)]:a.write_csv(OUT/(name+'.csv'),data)
    ok=all(x['qualified'] for x in dual+pairing);a.write(OUT/'summary.json',dict(candidate_promoted=ok,finite_qualified=sum(r['qualified']=='True' for r in dc),finite_count=16,full_curve_released=ok))
    d.matrix.freeze(OUT/'validation_evidence.json',[OUT/'freeze.json',OUT/'response_evidence.json',OUT/'summary.json']+list(OUT.glob('*.csv'))+[x for c in a.read(OUT/'contract.json')['cases'] for x in (LOCAL/c['key']).rglob('*') if x.is_file()])
    print('Finite mobility candidate promoted:',ok,flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('prepare','preflight','small','response','finite','analyze'));action=parser.parse_args().action
    if action in ('small','finite'):run(action=='finite')
    else:globals()[action]()
