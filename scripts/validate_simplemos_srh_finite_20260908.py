"""Matched local SRH volume candidates at qualified manual-constant states.

Four controls, individual and joint sources, independently calibrated before
finite substitution. All artifacts are new; original acceptance stays frozen.
"""
import argparse,copy,math
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from scipy.sparse import coo_matrix,diags
from scipy.sparse.linalg import splu
import build_simplemos_srh_finite_20260908 as b
import continue_simplemos_constants_finite_20260908 as continuation

c=b.c;a=c.a;d=c.d;v=c.v;p=c.p;LOCAL=b.LOCAL;OUT=b.OUT;RUNNER=b.RUNNER
AXES=('node_792','node_1057','joint')

def env(case,axis='zero',alpha=0.):
    e=c.env()
    if axis!='zero':e.update(VELA_CANDIDATE_SRH_VOLUMES=str(LOCAL/'volumes'/(case['device']+'_'+axis+'.txt')),VELA_CANDIDATE_SRH_ALPHA=format(alpha,'.17g'))
    return e
def execute(path,case,axis='zero',alpha=0.):return v.execute(path,RUNNER,env(case,axis,alpha))

def prepare():
    for f in (OUT/'build_evidence.json',c.OUT/'validation_evidence.json',continuation.OUT/'validation_evidence.json'):a.verify(f)
    assert a.read(continuation.OUT/'summary.json')['qualified']
    geometries={dev:c.s.prior.old.l.g.geometry(dev) for dev in ('n19','n23')};maps={};spatial=[]
    for dev in ('n19','n23'):
        geo,xy,el,vol,*_=geometries[dev];maps[dev]={}
        for tag in (792,1057):
            target=np.array(geometries['n23'][0].coords[tag]);dist=np.array([np.linalg.norm(np.array(geo.coords[i])-target) for i in range(geo.count)]);node=int(np.argmin(dist))
            assert node in geo.free and dist[node]<=1e-8 and np.sort(dist)[1]>1e-5
            ratio=float(vol['Si'][node]/geo.volumes['all_cell'][node]);assert ratio>0
            maps[dev][str(tag)]=dict(node=node,ratio=ratio)
            spatial.append(dict(device=dev,tag=tag,node=node,x_um=geo.coords[node][0],y_um=geo.coords[node][1],mapping_distance_um=float(dist[node]),second_nearest_distance_um=float(np.sort(dist)[1]),original_volume_m2=float(geo.volumes['all_cell'][node]),signed_Si_volume_m2=float(vol['Si'][node]),ratio=ratio))
        for axis,tags in [('node_792',[792]),('node_1057',[1057]),('joint',[792,1057])]:
            path=LOCAL/'volumes'/(dev+'_'+axis+'.txt');path.parent.mkdir(parents=True,exist_ok=True)
            assert not path.exists();path.write_text(f'{geo.count} {len(tags)}\n'+''.join(f"{maps[dev][str(tag)]['node']} {maps[dev][str(tag)]['ratio']:.17g}\n" for tag in tags))
    cs=c.cases();files=[Path(__file__).resolve(),OUT/'build_evidence.json',c.OUT/'validation_evidence.json',continuation.OUT/'validation_evidence.json',RUNNER]+list((LOCAL/'volumes').glob('*.txt'))
    for case in cs:
        case['mapped_nodes']=maps[case['device']];case['srh_probes']=[];case['srh_jobs']=[]
        cfg=c.old.cfg(case);assert cfg['solver']['recombination']==['srh'] and cfg['solver']['impact_ionization']['model']=='none'
        assert cfg['solver'].get('carrier_diagonal_floor',{}).get('scale',0.)==0.
        seed=c.LOCAL/case['key']/'manual_constants/state.csv';files.append(seed)
        for axis in ('zero',)+AXES:
            dest=LOCAL/case['key']/'fixed'/axis;paths=v.probes(cfg,dest,seed);selected=['functional','edges']
            if axis=='zero':selected+=['terms']
            for label,sim in [('jacobian','parameter_jacobian'),('jvp','joint_geometry_jvp')]:
                deck=copy.deepcopy(a.read(dest/'functional.json'));deck.update(simulation_type=sim,output_csv=str(dest/(label+'.csv')));a.write(dest/(label+'.json'),deck);paths.append(dest/(label+'.json'));selected.append(label)
            files+=paths
            case['srh_probes'] += [dict(path=str(x),axis=axis,alpha=0. if axis=='zero' else 1.) for x in paths if x.stem in selected]
            jobs=[('zero',0.,seed)] if axis=='zero' else [(name,alpha,seed) for name,alpha in [('plus_full',.001),('minus_full',-.001),('plus_half',.0005),('minus_half',-.0005),('finite',1.)]]
            if axis!='zero':jobs.append(('finite_native_seed',1.,continuation.LOCAL/case['key']/'full/state.csv'));files.append(continuation.LOCAL/case['key']/'full/state.csv')
            for label,alpha,initial in jobs:
                dest=LOCAL/case['key']/axis/label;deck=copy.deepcopy(cfg);deck.update(state_file=str(initial),output_state_file=str(dest/'state.csv'))
                a.write(dest/'config.json',deck);v.post_config(deck,dest);files += [dest/(n+'.json') for n in ('config','all_row','acceptance_edges')]+v.probes(deck,dest/'post',dest/'state.csv')
                case['srh_jobs'].append(dict(axis=axis,label=label,alpha=alpha,dest=str(dest)))
    a.write_csv(OUT/'spatial_scope.csv',spatial)
    a.write(OUT/'contract.json',dict(cases=cs,small_DC=52,finite_DC=24,
        parameter='V_i(alpha)=V_original_i*[1+alpha*(V_signed_Si_i/V_original_i-1)] at two spatially matched nodes only; individual and joint arms.',
        constants=c.NATIVE,source='Geometry ratios from exported/independently reconstructed native Si box measures; no field or current fitting.',
        scope='All four Vg=.8 controls; n23 nodes 792/1057 and n19 spatial matches 791/1056. Node match <=1e-8 um and next candidate farther than 1e-5 um.',
        preflight_gates=dict(zero_bit_identity=True,source_relative=1e-5,linear_relative=1e-8,weak_J_scaling_relative=1e-12,strong_JVP_relative=1e-4),
        calibration='Full physical potential response vector (psi,phin,phip) on all 907 free Si nodes, with independently solved full J tangent. Local fields reported separately.',
        response_gates=dict(prediction_relative=1e-3,two_amplitude_relative=1e-3,even_over_odd=.01,signal_over_zero_drift=100,zero_Id_drift_dex=1e-5),
        DC_gates=a.read(c.old.OUT/'contract.json')['unchanged_DC_gates'],
        finite='Same DC/dual/port gates; full field errors and selected-node changes reported without declaring global physical promotion or replacing Id-error gates.',
        native_response='Native SRH same-source self-consistent calibration not inferred from this Vela experiment.',defaults_changed=False))
    d.matrix.freeze(OUT/'freeze.json',files+[OUT/'contract.json',OUT/'spatial_scope.csv']);print('Frozen SRH 52 small +24 finite DC',flush=True)

def matrix(path,N):
    rows=a.rows(path);return coo_matrix(([float(x['value']) for x in rows],([int(x['row']) for x in rows],[int(x['column']) for x in rows])),shape=(3*N,3*N)).tocsc()

def preflight():
    a.verify(OUT/'freeze.json');checks=[];jvps=[];derivatives=[]
    for case in a.read(OUT/'contract.json')['cases']:
        for job in case['srh_probes']:
            status=execute(Path(job['path']),case,job['axis'],job['alpha']);assert status['exit_code']==0,(job,status)
        geo,mask=v.m.previous.prior.support(case);N=geo.count;ids=np.flatnonzero(mask);root=LOCAL/case['key']/'fixed';zero=root/'zero'
        R0=c.old.R(zero,N);previous=c.old.R(c.LOCAL/case['key']/'manual_constants/post',N);terms=d.ordered(zero/'terms.csv',N)
        J=matrix(zero/'jacobian.csv',N);scale=1./np.maximum(np.asarray(abs(J).max(axis=1).toarray()).ravel(),1e-300);lu=splu((diags(scale)@J).tocsc());wide=J.tocsr().astype(np.longdouble)
        potential=a.read(zero/'jacobian.status.json')['potential_scale_V'];baseedges=a.rows(zero/'edges.csv')
        r,_=c.s.prior.old.jc.jvp_metrics(zero/'jvp.csv',case['key'],'zero');jvps+=r
        assert all(x['qualified'] for x in r)
        for axis in AXES:
            selected=['792','1057'] if axis=='joint' else [axis[5:]];ratios={case['mapped_nodes'][tag]['node']:case['mapped_nodes'][tag]['ratio'] for tag in selected};source=np.zeros(3*N)
            for node,ratio in ratios.items():
                for block,car in ((1,'electron'),(2,'hole')):source[block*N+node]=float(terms[node][car+'_recombination'])*(ratio-1.)
            actual=c.old.R(root/axis,N)-R0;source_error=np.linalg.norm(actual-source)/max(np.linalg.norm(source),1e-300)
            changedJ=matrix(root/axis/'jacobian.csv',N);weak=0.
            for out,inp in ((1,2),(2,1)):
                base=J[out*N:(out+1)*N,inp*N:(inp+1)*N].toarray();new=changedJ[out*N:(out+1)*N,inp*N:(inp+1)*N].toarray();expected=base.copy()
                for node,ratio in ratios.items():expected[node,:]*=ratio
                nonzero=abs(expected)>0
                weak=max(weak,float(np.max(abs(new[nonzero]/expected[nonzero]-1))) if np.any(nonzero) else 0.)
                assert np.all(new[~nonzero]==0.)
            dx=lu.solve(-scale*source)
            for _ in range(4):dx-=lu.solve(scale*np.asarray(wide@dx.astype(np.longdouble)+source.astype(np.longdouble),dtype=float))
            linear=float(np.linalg.norm(np.asarray(wide@dx.astype(np.longdouble)+source.astype(np.longdouble),dtype=float))/max(np.linalg.norm(source),1e-300))
            for i in ids:derivatives.append(dict(key=case['key'],axis=axis,node_id=int(i),psi_V=potential*dx[i],phin_V=potential*dx[N+i],phip_V=potential*dx[2*N+i]))
            jr,_=c.s.prior.old.jc.jvp_metrics(root/axis/'jvp.csv',case['key'],axis);jvps+=jr
            edges=a.rows(root/axis/'edges.csv');unchanged=all(x[k]==y[k] for x,y in zip(edges,baseedges) for k in ('electron_flux','hole_flux','electron_mobility_m2_V_s','hole_mobility_m2_V_s','couple_m'))
            st=a.read(root/axis/'functional.status.json');port=abs(st['current_A_per_um']/st['contact_current_extractor_A_per_um']-1)
            passed=np.array_equal(R0,previous) and source_error<=1e-5 and linear<=1e-8 and weak<=1e-12 and unchanged and port<=1e-8 and all(x['qualified'] for x in jr)
            checks.append(dict(key=case['key'],axis=axis,zero_bit_identity=bool(np.array_equal(R0,previous)),source_relative=float(source_error),linear_relative=linear,weak_J_scaling_relative=weak,transport_bit_identity=unchanged,port_relative=port,qualified=bool(passed)))
        print('SRH preflight',case['key'],[x['qualified'] for x in checks if x['key']==case['key']],flush=True)
    for name,data in [('preflight',checks),('jvp',jvps),('state_derivative',derivatives)]:a.write_csv(OUT/(name+'.csv'),data)
    d.matrix.freeze(OUT/'preflight_evidence.json',[OUT/'freeze.json']+[OUT/(n+'.csv') for n in ('preflight','jvp','state_derivative')]+[x for case in c.cases() for x in (LOCAL/case['key']/'fixed').rglob('*') if x.is_file()])

def run(finite=False):
    a.verify(OUT/'freeze.json');a.verify(OUT/'preflight_evidence.json');assert all(x['qualified']=='True' for x in a.rows(OUT/'preflight.csv'))
    if finite:
        a.verify(OUT/'response_evidence.json');assert all(x['qualified']=='True' for x in a.rows(OUT/'response.csv'))
    def one(case):
        rows=[]
        for job in case['srh_jobs']:
            if job['label'].startswith('finite')!=finite:continue
            dest=Path(job['dest']);status=execute(dest/'config.json',case,job['axis'],job['alpha']);row=dict(key=case['key'],axis=job['axis'],label=job['label'],qualified=False,failure=status.get('failure_reason',''))
            if (dest/'state.csv').exists():
                for n in ('all_row','acceptance_edges'):execute(dest/(n+'.json'),case,job['axis'],job['alpha'])
                row.update(c.s.prior.old.qualify(case,dest))
                if finite:
                    for n in ('functional','edges','terms'):execute(dest/'post'/(n+'.json'),case,job['axis'],job['alpha'])
            rows.append(row);print(case['key'],job['axis'],job['label'],row['qualified'],row.get('max_row_ratio'),flush=True)
        return rows
    with ThreadPoolExecutor(max_workers=2) as pool:rows=[r for group in pool.map(one,a.read(OUT/'contract.json')['cases']) for r in group]
    keys=list(dict.fromkeys(k for r in rows for k in r));a.write_csv(OUT/('finite_dc.csv' if finite else 'small_dc.csv'),[{k:r.get(k,'') for k in keys} for r in rows])

def physical_delta(state,reference,ids):
    return np.array([[float(v.m.previous.prior.physical(state[i],field)-v.m.previous.prior.physical(reference[i],field)) for i in ids] for field in ('psi','phin','phip')])

def response():
    dc=a.rows(OUT/'small_dc.csv');derivs=a.rows(OUT/'state_derivative.csv');results=[]
    for case in a.read(OUT/'contract.json')['cases']:
        geo,mask=v.m.previous.prior.support(case);ids=np.flatnonzero(mask);root=LOCAL/case['key'];original=d.ordered(c.LOCAL/case['key']/'manual_constants/state.csv',geo.count)
        zero=d.ordered(root/'zero/zero/state.csv',geo.count);noise=float(np.linalg.norm(physical_delta(zero,original,ids)))
        rows={(r['axis'],r['label']):r for r in dc if r['key']==case['key']};zeroid=float(rows[('zero','zero')]['current_A_per_um']);oldid=a.read(c.LOCAL/case['key']/'manual_constants/config.status.json')['contact_currents_A_per_um']['drain'];dex=abs(math.log10(zeroid/oldid))
        for axis in AXES:
            dr=[r for r in derivs if r['key']==case['key'] and r['axis']==axis];assert [int(r['node_id']) for r in dr]==list(ids)
            tangent=np.array([[float(r[f+'_V']) for r in dr] for f in ('psi','phin','phip')]);delta={}
            for label in ('plus_full','minus_full','plus_half','minus_half'):delta[label]=physical_delta(d.ordered(root/axis/label/'state.csv',geo.count),zero,ids)
            odd={k:(delta['plus_'+k]-delta['minus_'+k])/2 for k in ('full','half')};lin=float(np.linalg.norm(odd['full']-2*odd['half'])/max(np.linalg.norm(odd['full']),1e-300))
            for amp,alpha in [('full',.001),('half',.0005)]:
                signal=float(np.linalg.norm(odd[amp]));error=float(np.linalg.norm(odd[amp]-alpha*tangent)/max(np.linalg.norm(alpha*tangent),1e-300));even=float(np.linalg.norm((delta['plus_'+amp]+delta['minus_'+amp])/2)/max(signal,1e-300));snr=signal/max(noise,1e-300)
                passed=all(rows[(axis,sgn+'_'+amp)]['qualified']=='True' for sgn in ('plus','minus')) and rows[('zero','zero')]['qualified']=='True' and error<=1e-3 and lin<=1e-3 and even<=.01 and snr>=100 and dex<=1e-5
                results.append(dict(key=case['key'],axis=axis,amplitude=amp,prediction_relative=error,two_amplitude_relative=lin,even_over_odd=even,signal_over_zero_drift=snr,zero_Id_drift_dex=dex,odd_norm_V=signal,qualified=passed))
    a.write_csv(OUT/'response.csv',results);d.matrix.freeze(OUT/'response_evidence.json',[OUT/'preflight_evidence.json',OUT/'small_dc.csv',OUT/'response.csv']+[x for case in a.read(OUT/'contract.json')['cases'] for job in case['srh_jobs'] if not job['label'].startswith('finite') for x in Path(job['dest']).rglob('*') if x.is_file()]);print('Qualified small SRH',sum(r['qualified'] for r in results),'/',len(results),flush=True)

def analyze():
    dc=a.rows(OUT/'finite_dc.csv');duals=[];fields=[];local=[];ports=[];currents=[]
    for case in a.read(OUT/'contract.json')['cases']:
        geo,mask=v.m.previous.prior.support(case);baseid=a.read(c.LOCAL/case['key']/'manual_constants/config.status.json')['contact_currents_A_per_um']['drain'];selected=[case['mapped_nodes'][tag]['node'] for tag in ('792','1057')]+[1000,1009]
        for axis in AXES:
            tags=['792','1057'] if axis=='joint' else [axis[5:]];ratios={case['mapped_nodes'][tag]['node']:case['mapped_nodes'][tag]['ratio'] for tag in tags};states=[];ids=[]
            for label in ('finite','finite_native_seed'):
                dest=LOCAL/case['key']/axis/label;row=next(r for r in dc if r['key']==case['key'] and r['axis']==axis and r['label']==label);Id=float(row['current_A_per_um']);ids.append(Id);states.append(d.ordered(dest/'state.csv',geo.count))
                met,loc=c.field_metrics(case,axis+'/'+label,dest,ratios,selected);fields+=met;local+=loc
                status=a.read(dest/'post/functional.status.json');error=max(abs(Id/status['current_A_per_um']-1),abs(status['current_A_per_um']/status['contact_current_extractor_A_per_um']-1));ports.append(dict(key=case['key'],axis=axis,label=label,relative=error,qualified=error<=1e-8))
                currents.append(dict(key=case['key'],device=case['device'],vd=case['vd'],axis=axis,label=label,Id_A_per_um=Id,native_relative_error=Id/case['native_Id_A_per_um']-1,delta_over_constant_baseline=Id/baseid-1,qualified=row['qualified']=='True'))
            diff=v.m.previous.prior.delta_states(*states,mask);ie=abs(ids[0]/ids[1]-1)
            duals.append(dict(key=case['key'],axis=axis,**diff,Id_relative=ie,qualified=max(diff[k] for k in ('psi_max_V','phin_max_V','phip_max_V'))<=1e-6 and diff['density_max_relative']<=1e-4 and ie<=1e-6))
    for name,data in [('dual',duals),('fields',fields),('local_fields',local),('ports',ports),('currents',currents)]:a.write_csv(OUT/(name+'.csv'),data)
    summary=dict(DC=len(dc),qualified_DC=sum(r['qualified']=='True' for r in dc),qualified_dual=sum(r['qualified'] for r in duals),qualified_ports=sum(r['qualified'] for r in ports),max_row_ratio=max(float(r['max_row_ratio']) for r in dc),qualified=all(r['qualified']=='True' for r in dc) and all(r['qualified'] for r in duals+ports),native_source_calibrated=False)
    a.write(OUT/'summary.json',summary);d.matrix.freeze(OUT/'validation_evidence.json',[OUT/'freeze.json',OUT/'preflight_evidence.json',OUT/'response_evidence.json',OUT/'summary.json']+list(OUT.glob('*.csv'))+[x for case in a.read(OUT/'contract.json')['cases'] for job in case['srh_jobs'] if job['label'].startswith('finite') for x in Path(job['dest']).rglob('*') if x.is_file()]);print(summary,flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('prepare','preflight','small','response','finite','analyze'));action=parser.parse_args().action
    if action in ('small','finite'):run(action=='finite')
    else:globals()[action]()
