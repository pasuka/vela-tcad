"""Same-parameter residual/tangent versus independent nonlinear DC calibration.

Constants are perturbed coherently in every compiled core translation unit.
SRH tests use a local field observable fixed before execution: a vanishing drain
response cannot qualify a minority-field source. No native simulator constants
are fitted, and no finite physical replacement is authorized by this driver.
"""
import argparse,copy,math
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from scipy.sparse import coo_matrix,diags
from scipy.sparse.linalg import splu
import build_simplemos_constants_srh_calibration_20260907 as b

s=b.s;a=b.a;d=b.d;v=b.v;p=b.p;LOCAL=b.LOCAL;OUT=b.OUT;RUNNER=b.RUNNER
AXES=('q_fixed_Vt','eps0','thermal_voltage')

def env(axis,alpha):
    e=v.environment()
    if axis=='q_fixed_Vt':e.update(VELA_CANDIDATE_Q_RELATIVE=str(alpha),VELA_CANDIDATE_KB_RELATIVE=str(alpha))
    elif axis=='eps0':e['VELA_CANDIDATE_EPS0_RELATIVE']=str(alpha)
    elif axis=='thermal_voltage':e['VELA_CANDIDATE_KB_RELATIVE']=str(alpha)
    elif axis.startswith('srh_'):e.update(VELA_CANDIDATE_SRH_NODE=axis[4:],VELA_CANDIDATE_SRH_RELATIVE=str(alpha))
    elif axis!='zero':raise ValueError(axis)
    return e

def execute(path,axis='zero',alpha=0.):return v.execute(path,RUNNER,env(axis,alpha))
def cfg(c):return s.config(c,'generated')
def cases():return [c for c in s.prior.cases() if c['vg']==.8]

def prepare():
    a.verify(OUT/'build_evidence.json');a.verify(s.OUT/'validation_evidence.json')
    assert a.read(s.OUT/'summary.json')['qualified_DC']==32
    files=[Path(__file__).resolve(),OUT/'build_evidence.json',s.OUT/'validation_evidence.json',RUNNER]
    cs=cases()
    for c in cs:
        state=s.LOCAL/c['key']/'generated/state.csv';c['seed']=str(state);c['calibration_axes']=list(AXES)
        if c['device']=='n23':c['calibration_axes']+=['srh_792','srh_1057']
        c['fixed_jobs']=[];c['dc_jobs']=[]
        for axis in ['zero']+c['calibration_axes']:
            amps=[('zero',0.)] if axis=='zero' else [(label,alpha) for label,alpha in [('p1',1e-5),('m1',-1e-5),('p2',5e-6),('m2',-5e-6)]]
            if axis.startswith('srh_'):amps=[('p1',.001),('m1',-.001)]
            for label,alpha in amps:
                dest=LOCAL/c['key']/'fixed'/axis/label;paths=v.probes(cfg(c),dest,state)
                selected=['functional']
                if axis=='zero':
                    selected+=['terms','adjoint']
                    j=copy.deepcopy(a.read(dest/'functional.json'));j.update(simulation_type='parameter_jacobian',output_csv=str(dest/'jacobian.csv'));a.write(dest/'jacobian.json',j);paths.append(dest/'jacobian.json');selected.append('jacobian')
                files+=paths
                c['fixed_jobs'] += [dict(path=str(x),axis=axis,alpha=alpha) for x in paths if x.stem in selected]
            amp=.001 if axis.startswith('srh_') else .0001
            amps=[('zero',0.)] if axis=='zero' else [('plus_full',amp),('minus_full',-amp),('plus_half',amp/2),('minus_half',-amp/2)]
            for label,alpha in amps:
                dest=LOCAL/c['key']/'dc'/axis/label;deck=cfg(c);deck.update(state_file=str(state),output_state_file=str(dest/'state.csv'))
                a.write(dest/'config.json',deck);v.post_config(deck,dest)
                files += [dest/(name+'.json') for name in ('config','all_row','acceptance_edges')]
                c['dc_jobs'].append(dict(axis=axis,label=label,alpha=alpha,dest=str(dest)))
        files.append(state)
    a.write(OUT/'contract.json',dict(cases=cs,DC=sum(len(c['dc_jobs']) for c in cs),
        scope='Four Vg=.8 controls, two NWell and two Vd; node 792/1057 SRH local fields at both n23 Vd. Masetti/SRH-only profile unchanged.',
        directions={'q_fixed_Vt':'q and kb both multiply by 1+alpha, keeping Vt fixed; every q occurrence changes coherently.',
                    'eps0':'eps0 multiplies by 1+alpha in all geometry/scaling/material paths.',
                    'thermal_voltage':'kb multiplies by 1+alpha at fixed q,T,ni and material inputs.',
                    'srh':'At the selected node only, SRH source volume and its state derivatives multiply by 1+alpha; independent Poisson charge volumes fixed. Avalanche and diagonal floors off.'},
        constants_amplitudes=[.0001,.00005],srh_amplitudes=[.001,.0005],
        prediction='Fixed-physical-state residual parameter derivative plus explicit port derivative, checked at two steps. Full base J tangent and calibrated terminal adjoint. SRH source derivative computed directly from independent solved source terms.',
        observables={'constants':'drain current','srh_792':'phin at node 792','srh_1057':'phip at node 1057'},
        gates=dict(prediction_relative=1e-3,two_amplitude_relative=1e-3,even_over_odd=.01,signal_over_zero_drift=100,zero_Id_drift_dex=1e-5,linear_residual_relative=1e-8),
        unchanged_DC_gates=a.read(s.OUT/'contract.json')['gates'],native_new_runs=False,
        boundary='Vela same-parameter self-consistent calibration only. No claim that this changes or independently calibrates native Sentaurus global constants or SRH source definitions. No finite substitution in this driver.'))
    d.matrix.freeze(OUT/'freeze.json',files+[OUT/'contract.json'])
    print('Frozen',sum(len(c['dc_jobs']) for c in cs),'DC',flush=True)

def R(path,n):return d.array(d.ordered(path/'residual.csv',n),('psi_residual','phin_residual','phip_residual')).reshape(-1)
def field(state,node,kind):return float(v.m.previous.prior.physical(state[node],kind))

def preflight():
    a.verify(OUT/'freeze.json');records=[];preds=[]
    for c in a.read(OUT/'contract.json')['cases']:
        for job in c['fixed_jobs']:
            status=execute(Path(job['path']),job['axis'],job['alpha']);assert status['exit_code']==0,(job,status)
        root=LOCAL/c['key']/'fixed';base=root/'zero/zero';geo,mask=v.m.previous.prior.support(c);N=geo.count
        zero=R(base,N);old=R(s.LOCAL/c['key']/'generated/post',N)
        status=a.read(base/'functional.status.json');previous=a.read(s.LOCAL/c['key']/'generated/post/functional.status.json')
        # Runtime zero must retain the production operator and port at the same state.
        zero_delta=np.linalg.norm(zero-old)/max(np.linalg.norm(old),1e-300)
        zero_Id=abs(status['current_A_per_um']/previous['current_A_per_um']-1)
        lam=d.array(d.ordered(base/'adjoint.csv',N),('lambda_poisson','lambda_electron','lambda_hole')).reshape(-1)
        rows=a.rows(base/'jacobian.csv');J=coo_matrix(([float(x['value']) for x in rows],([int(x['row']) for x in rows],[int(x['column']) for x in rows])),shape=(3*N,3*N)).tocsc()
        scaling=1./np.maximum(np.asarray(abs(J).max(axis=1).toarray()).ravel(),1e-300)
        lu=splu((diags(scaling)@J).tocsc());Jwide=J.tocsr().astype(np.longdouble)
        potential=a.read(base/'jacobian.status.json')['potential_scale_V']
        terms=d.ordered(base/'terms.csv',N)
        for axis in c['calibration_axes']:
            ar=root/axis;direct=0.;source_check=0.
            if axis.startswith('srh_'):
                node=int(axis[4:]);source=np.zeros(3*N)
                for block,carrier in ((1,'electron'),(2,'hole')):source[block*N+node]=float(terms[node][carrier+'_recombination'])
                actual=(R(ar/'p1',N)-R(ar/'m1',N))/.002
                source_check=np.linalg.norm(actual-source)/max(np.linalg.norm(source),1e-300)
                observable='phin' if node==792 else 'phip';block=1 if node==792 else 2
            else:
                derivatives=[(R(ar/f'p{i}',N)-R(ar/f'm{i}',N))/(2*h) for i,h in ((1,1e-5),(2,5e-6))]
                source=derivatives[1];source_check=np.linalg.norm(derivatives[0]-source)/max(np.linalg.norm(source),1e-300)
                currents=[(a.read(ar/f'p{i}/functional.status.json')['current_A_per_um']-a.read(ar/f'm{i}/functional.status.json')['current_A_per_um'])/(2*h) for i,h in ((1,1e-5),(2,5e-6))]
                direct=currents[1];observable='Id';node=-1;block=-1
            dx=lu.solve(-scaling*source)
            for _ in range(4):
                residual=np.asarray(Jwide@dx.astype(np.longdouble)+source.astype(np.longdouble),dtype=float);dx-=lu.solve(scaling*residual)
            residual=np.asarray(Jwide@dx.astype(np.longdouble)+source.astype(np.longdouble),dtype=float)
            lin=np.linalg.norm(residual)/max(np.linalg.norm(source),1e-300)
            current=direct-math.fsum(float(x*y) for x,y in zip(lam,source))
            prediction=current if observable=='Id' else potential*dx[block*N+node]
            preds.append(dict(key=c['key'],axis=axis,observable=observable,node=node,unit_prediction=prediction,unit_Id_prediction_A_per_um=current,unit_direct_Id_A_per_um=direct))
            records.append(dict(key=c['key'],axis=axis,zero_residual_bit_identity=bool(np.array_equal(zero,old)),zero_residual_relative=zero_delta,zero_Id_relative=zero_Id,parameter_source_relative=source_check,linear_residual_relative=lin,
                                qualified=np.array_equal(zero,old) and zero_Id<=1e-8 and source_check<=1e-5 and lin<=1e-8))
        print('Parameter preflight',c['key'],[(x['axis'],x['qualified']) for x in records if x['key']==c['key']],flush=True)
    a.write_csv(OUT/'preflight.csv',records);a.write_csv(OUT/'predictions.csv',preds)
    d.matrix.freeze(OUT/'preflight_evidence.json',[OUT/'freeze.json',OUT/'preflight.csv',OUT/'predictions.csv']+[x for c in cases() for x in (LOCAL/c['key']/'fixed').rglob('*') if x.is_file()])

def run():
    a.verify(OUT/'freeze.json');a.verify(OUT/'preflight_evidence.json')
    assert all(x['qualified']=='True' for x in a.rows(OUT/'preflight.csv'))
    def one(c):
        rows=[]
        for j in c['dc_jobs']:
            dest=Path(j['dest']);status=execute(dest/'config.json',j['axis'],j['alpha'])
            r=dict(key=c['key'],axis=j['axis'],label=j['label'],alpha=j['alpha'],qualified=False,failure=status.get('failure_reason',''))
            if (dest/'state.csv').exists():
                for name in ('all_row','acceptance_edges'):execute(dest/(name+'.json'),j['axis'],j['alpha'])
                r.update(s.prior.old.qualify(c,dest))
            rows.append(r);print(c['key'],j['axis'],j['label'],r['qualified'],r.get('max_row_ratio'),flush=True)
        return rows
    with ThreadPoolExecutor(max_workers=2) as pool:data=[r for rows in pool.map(one,a.read(OUT/'contract.json')['cases']) for r in rows]
    keys=list(dict.fromkeys(k for r in data for k in r));a.write_csv(OUT/'dc.csv',[{k:r.get(k,'') for k in keys} for r in data])

def analyze():
    a.verify(OUT/'freeze.json');dc=a.rows(OUT/'dc.csv');pred=a.rows(OUT/'predictions.csv');response=[];local=[]
    for c in a.read(OUT/'contract.json')['cases']:
        geo,mask=v.m.previous.prior.support(c);root=LOCAL/c['key']/'dc';z=next(x for x in dc if x['key']==c['key'] and x['axis']=='zero')
        states={('zero','zero'):d.ordered(root/'zero/zero/state.csv',geo.count)}
        original=d.ordered(Path(c['seed']),geo.count);oldId=a.read(s.LOCAL/c['key']/'generated/config.status.json')['contact_currents_A_per_um']['drain']
        zeroId=float(z['current_A_per_um']);zeroDex=abs(math.log10(zeroId/oldId))
        for pr in [r for r in pred if r['key']==c['key']]:
            axis=pr['axis'];node=int(pr['node']);obs=pr['observable'];prediction=float(pr['unit_prediction'])
            rows={r['label']:r for r in dc if r['key']==c['key'] and r['axis']==axis};values={}
            for label,r in rows.items():
                state=d.ordered(root/axis/label/'state.csv',geo.count);states[(axis,label)]=state
                values[label]=float(r['current_A_per_um']) if obs=='Id' else field(state,node,obs)
                if node>=0:
                    terms=d.ordered(root/axis/label/'all_row.csv',geo.count)
                    local.append(dict(key=c['key'],axis=axis,label=label,node=node,alpha=r['alpha'],psi_V=field(state,node,'psi'),phin_V=field(state,node,'phin'),phip_V=field(state,node,'phip'),electrons_m3=state[node]['electrons_m3'],holes_m3=state[node]['holes_m3'],electron_SRH_source=terms[node]['electron_recombination'],hole_SRH_source=terms[node]['hole_recombination'],Id_A_per_um=r['current_A_per_um']))
            center=zeroId if obs=='Id' else field(states[('zero','zero')],node,obs)
            initial=oldId if obs=='Id' else field(original,node,obs);drift=abs(center-initial)
            odd={k:(values['plus_'+k]-values['minus_'+k])/2 for k in ('full','half')}
            linearity=abs(odd['full']/(2*odd['half'])-1) if odd['half'] else math.inf
            for amp in ('full','half'):
                alpha=float(rows['plus_'+amp]['alpha']);error=abs(odd[amp]/(alpha*prediction)-1) if prediction else math.inf
                even=abs((values['plus_'+amp]+values['minus_'+amp])/2-center)/abs(odd[amp]) if odd[amp] else math.inf
                snr=abs(odd[amp])/max(drift,1e-300);signs=(values['plus_'+amp]-center)*prediction>0 and (values['minus_'+amp]-center)*prediction<0
                qualified=z['qualified']=='True' and all(rows[k+'_'+amp]['qualified']=='True' for k in ('plus','minus')) and error<=1e-3 and linearity<=1e-3 and even<=.01 and snr>=100 and zeroDex<=1e-5 and signs
                response.append(dict(key=c['key'],axis=axis,observable=obs,node=node,amplitude=amp,alpha=alpha,prediction=alpha*prediction,actual_odd=odd[amp],prediction_relative_error=error,two_amplitude_relative=linearity,even_over_odd=even,signal_over_zero_drift=snr,zero_Id_drift_dex=zeroDex,qualified=qualified))
    a.write_csv(OUT/'response.csv',response);a.write_csv(OUT/'local_fields.csv',local)
    summary=dict(DC=len(dc),qualified_DC=sum(r['qualified']=='True' for r in dc),response_checks=len(response),qualified_response=sum(r['qualified'] for r in response),native_same_perturbation_done=False,finite_replacement_done=False)
    a.write(OUT/'summary.json',summary)
    d.matrix.freeze(OUT/'validation_evidence.json',[OUT/'freeze.json',OUT/'preflight_evidence.json',OUT/'summary.json']+list(OUT.glob('*.csv'))+[x for c in cases() for x in (LOCAL/c['key']/'dc').rglob('*') if x.is_file()]);print(summary,flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('prepare','preflight','run','analyze'));globals()[parser.parse_args().action]()
