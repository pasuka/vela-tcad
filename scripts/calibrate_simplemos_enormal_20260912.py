"""Native Enormal observer and unfitted spatial formation calibration.

The old inferred PhuMob floors are diagnostic controls only. No fitting here,
and no promotion of those floors to production or relaxation of old gates.
"""
import argparse, math, tarfile
from pathlib import Path
from collections import defaultdict
from functools import lru_cache
from decimal import localcontext
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from scipy.spatial import cKDTree
import simplemos_phumob_curves_20260912 as c
import audit_simplemos_masetti_native_geometry_20260907 as g
import audit_simplemos_phumob_cross_derivatives_20260908 as h
import localize_simplemos_phumob_screening_floor_20260909 as f

R=Path(__file__).resolve().parents[1];a,d=c.a,c.d
L=R/'build-release/enormal_restore_20260912'
O=R/'reference_tcad/simplemos_sentaurus2022/enormal_restore_20260912'
FLOORS=f.OUT/'floor_inference.csv'
PARAMS=((4.75e7,580.,.125,5.82e14,5.82e30),(9.925e6,2947.,.0317,2.0546e14,2.0546e30))

def unpack(expected):
    assert a.sha(L/'results.tgz')==expected
    raw=L/'native_raw';assert not raw.exists()
    with tarfile.open(L/'results.tgz') as t:
        for m in t.getmembers():
            assert (raw/m.name).resolve().is_relative_to(raw.resolve()) and not m.issym() and not m.islnk()
        t.extractall(raw,filter='data')
    for p in (L/'bundle').rglob('*'):
        if p.is_file():assert a.sha(p)==a.sha(raw/'bundle'/p.relative_to(L/'bundle'))
    assert a.sha(L/'pmi_vela_enormal_observer.C')==a.sha(raw/'pmi/pmi_vela_enormal_observer.C')
    points=[];exports=[]
    for j in a.read(O/'native_contract.json')['jobs']:
        root=raw/'bundle'/j['name'];code=int((root/'exit_code.txt').read_text());log=(root/'console.log').read_text(errors='replace')
        vals=c.n.exporter.pltrows(root/'native_des.plt');assert len(vals)==1
        v=vals[0];currents=[v[k+' TotalCurrent'] for k in ('drain','source','gate','substrate')]
        kcl=abs(math.fsum(currents))/max(abs(currents[0]),1e-300)
        bias=max(abs(v['gate OuterVoltage']-j['vg']),abs(v['drain OuterVoltage']-j['vd']))
        good=code==0 and 'T-2022.03-SP2' in log and 'Good Bye' in log and bias<=1e-10 and kcl<=1e-8
        points.append(dict(**j,exit_code=code,Id_A_per_um=currents[0],kcl_over_Id=kcl,bias_error_V=bias,qualified=good));assert good,points[-1]
        exports.append(dict(case=j['case'],index=j['index'],tdr=str(root/'final_des.tdr'),export=str(L/'exports'/j['name'])))
    a.write_csv(O/'native_points.csv',points);a.write(O/'export_contract.json',dict(jobs=exports))
    d.matrix.freeze(O/'native_evidence.json',[O/'native_freeze.json',O/'native_points.csv',O/'export_contract.json',L/'results.tgz']+[p for p in raw.rglob('*') if p.is_file()])

def export():
    a.verify(O/'native_evidence.json')
    with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(c.n.exporter.export_one,a.read(O/'export_contract.json')['jobs']))
    d.matrix.freeze(O/'export_evidence.json',[O/'native_evidence.json']+[p for p in (L/'exports').rglob('*') if p.is_file()])

def scalar(src,field,index='node_id'):
    return {int(r[index]):float(r['component0']) for r in a.rows(src/'fields'/field)}

@lru_cache(None)
def geometry(device):
    geo,xy,cells,vol,K,parts,info,mod=g.geometry(device)
    measures=g.box.parse_debug_block((g.LOCAL/'native_raw/bundle'/device/'MeasureCoefficients.debug').read_text(),'Measure');vp=info['measure_permutation']
    adj=defaultdict(list)
    for cid,cell in cells.items():
        for k in range(3):adj[tuple(sorted((cell['nodes'][k],cell['nodes'][(k+1)%3])))].append(cid)
    interfaces=[e for e,cs in adj.items() if 'Si' in {cells[k]['material'] for k in cs} and len({cells[k]['material'] for k in cs})>1]
    si={cid:cell for cid,cell in cells.items() if cell['material']=='Si'};nodes=sorted({k for cell in si.values() for k in cell['nodes']})
    distance={}
    for k in nodes:
        point=np.array(xy[k]);ds=[]
        for i,j in interfaces:
            p=np.array(xy[i]);v=np.array(xy[j])-p;t=np.clip(np.dot(point-p,v)/np.dot(v,v),0,1)
            ds.append(float(np.linalg.norm(point-p-t*v)))
        distance[k]=min(ds)
    stencil={};gd={};weights={};nodecells=defaultdict(list)
    for cid,cell in si.items():
        ns=cell['nodes'];pts=np.array([xy[k] for k in ns]);u=pts[1]-pts[0];v=pts[2]-pts[0];det=u[0]*v[1]-v[0]*u[1]
        # Differences first preserve constant-potential and constant-distance invariance.
        stencil[cid]=np.array([[v[1],-v[0]],[-u[1],u[0]]])/det
        gd[cid]=np.array([distance[ns[1]]-distance[ns[0]],distance[ns[2]]-distance[ns[0]]])@stencil[cid]
        m=[measures[cid]['values'][vp[k]] for k in range(3)];weights[cid]=[w/math.fsum(m) for w in m]
        for k in ns:nodecells[k].append(cid)
    return dict(xy=xy,cells=si,distance=distance,gradient=stencil,gd=gd,weights=weights,nodecells=nodecells)

@lru_cache(None)
def minima():
    with localcontext() as ctx:
        ctx.prec=80;return [tuple(map(float,h.gminimum(i))) for i in range(2)]

def bulk(state,b,floor=None):
    nd,na,n,p=(state[k] for k in ('nd','na','n','p'))
    nd*=1+nd*nd/(.21*nd*nd+4e20**2);na*=1+na*na/(.5*na*na+7.2e20**2)
    mass=(1,1.258)[b];P=1/(2.459*(nd+na+(p if b==0 else n))**(2/3)/3.97e13+3.828*(n+p)/(1.36e20*mass))
    lo,gg=minima()[b];gg=gg if floor is None else floor
    if P>=lo:gg=max(gg,1-.89233*(.41372+P*(1/mass)**.28227)**(-.19778)+.005978*(P*mass**.72169)**(-1.80618))
    return f.mobility(state,('e','h')[b],gg)

def inverse_surface(F,distance_um,impurity,b):
    B,C,lamb,delta,eta=PARAMS[b]
    return math.exp(-distance_um/.01)*(F/(B+C*(impurity+1)**lamb*F**(2/3))+F*F/delta+F**3/eta)

def identity(base,obs):
    result=[]
    for p in (base/'fields').glob('*.csv'):
        br=a.rows(p);orr=a.rows(obs/'fields'/p.name);assert len(br)==len(orr)
        keys=[k for k in br[0] if k.startswith('component')];da=[];dr=[]
        for x,y in zip(br,orr):
            assert {k:v for k,v in x.items() if k not in keys}=={k:v for k,v in y.items() if k not in keys}
            for k in keys:
                u,w=float(x[k]),float(y[k]);da.append(abs(u-w));dr.append(abs(u-w)/max(abs(u),abs(w),1e-300))
        result.append(dict(field=p.name,count=len(br),max_absolute=max(da),max_relative=max(dr)))
    return result

def analyze():
    a.verify(O/'export_evidence.json');points=a.rows(O/'native_points.csv')
    floors={r['carrier']:float(r['inferred_effective_G_floor']) for r in a.rows(FLOORS)}
    identities=[];samples=[];formation=[];spatial=[];summary=[]
    for job in a.read(O/'native_contract.json')['jobs']:
        if job['arm']!='observed':continue
        key=job['name'].removesuffix('_observed');src=L/'exports'/job['name'];base=L/'exports'/(key+'_baseline');geom=geometry(job['device'])
        assert geom['xy']=={int(r['id']):(float(r['x_um']),float(r['y_um'])) for r in a.rows(src/'nodes.csv')}
        ids=sorted(geom['xy']);tree=cKDTree([geom['xy'][k] for k in ids])
        fields={k:scalar(src,v+'_region0.csv') for k,v in [('nd','DonorConcentration'),('na','AcceptorConcentration'),('n','eDensity'),('p','hDensity'),('psi','ElectrostaticPotential')]}
        ident=identity(base,src);identities.extend(dict(key=key,**r) for r in ident)
        idb=float(next(r for r in points if r['name']==key+'_baseline')['Id_A_per_um']);ido=float(next(r for r in points if r['name']==job['name'])['Id_A_per_um'])
        nativefield={int(r['cell_id']):np.array([float(r['component0']),float(r['component1'])]) for r in a.rows(src/'fields/ElectricField_region0_cells.csv')}
        F={};field_errors=[]
        for cid,cell in geom['cells'].items():
            ns=cell['nodes'];psi=fields['psi'];ef=-np.array([psi[ns[1]]-psi[ns[0]],psi[ns[2]]-psi[ns[0]]])@geom['gradient'][cid]*1e4
            F[cid]=abs(float(ef@geom['gd'][cid]));field_errors.append(float(max(abs(ef-nativefield[cid]))))
        observed_nodes=set();samp=[]
        for car in ('e','h'):
            for r in a.rows(L/'native_raw/bundle'/job['name']/f'enormal_observer_{car}_Silicon_1.csv'):
                dist,idx=tree.query([float(r['x_um']),float(r['y_um'])]);k=ids[int(idx)];assert dist<1e-12 and k in fields['n'];observed_nodes.add(k)
                error,cid=min((abs(F[cid]-float(r['enorm_native'])),cid) for cid in geom['nodecells'][k])
                rec=dict(key=key,carrier=car,node=k,closest_cell=cid,field_error_V_per_cm=error,potential_error_V=abs(float(r['potential_V'])-fields['psi'][k]),n_relative=abs(float(r['n_native'])/fields['n'][k]-1),p_relative=abs(float(r['p_native'])/fields['p'][k]-1),distance_error_um=abs(float(r['interface_distance_um'])-geom['distance'][k]),compute_distance_error_cm=abs(float(r['dist_native'])-geom['distance'][k]*1e-4),enormal_native=float(r['enorm_native']))
                samples.append(rec);samp.append(rec)
        missing=set(fields['n'])-observed_nodes
        for b,car in enumerate(('e','h')):
            native=scalar(src,car+'Mobility_region0_cells.csv','cell_id')
            for mode,floor in [('mathematical_floor',None),('prior_diagnostic_floor',floors[car])]:
                mu={k:bulk({q:fields[q][k] for q in ('nd','na','n','p')},b,floor) for k in fields['n']};errs=[];miss_bound=[]
                for cid,cell in geom['cells'].items():
                    values=[]
                    for k in cell['nodes']:
                        inv=inverse_surface(F[cid],geom['distance'][k],fields['nd'][k]+fields['na'][k],b)
                        values.append(1/(1/mu[k]+inv))
                        if k in missing:miss_bound.append(mu[k]*inv)
                    val=math.fsum(w*v for w,v in zip(geom['weights'][cid],values));err=val/native[cid]-1;errs.append(abs(err))
                    formation.append(dict(key=key,carrier=car,mode=mode,cell=cid,field_V_per_cm=F[cid],native_mu=native[cid],formula_mu=val,relative_error=err))
                summary.append(dict(key=key,carrier=car,mode=mode,cells=len(errs),max_relative=max(errs),p95_relative=float(np.percentile(errs,95)),below_original_cell_gate=max(errs)<=1e-7,missing_node_surface_relative_bound=max(miss_bound,default=0)))
        spatial.append(dict(key=key,current_identity_relative=abs(ido/idb-1),field_identity_max_abs=max(r['max_absolute'] for r in ident),field_identity_max_relative=max(r['max_relative'] for r in ident),identity_fields=len(ident),observed_nodes=len(observed_nodes),silicon_nodes=len(fields['n']),missing_nodes=len(missing),cell_gradient_max_V_per_cm=max(field_errors),observer_field_max_V_per_cm=max(r['field_error_V_per_cm'] for r in samp),observer_potential_max_V=max(r['potential_error_V'] for r in samp),observer_n_max_relative=max(r['n_relative'] for r in samp),observer_p_max_relative=max(r['p_relative'] for r in samp),distance_max_um=max(r['distance_error_um'] for r in samp)))
    for name,rows in [('field_identity.csv',identities),('observer_samples.csv',samples),('cell_formation.csv',formation),('spatial_summary.csv',spatial),('formation_summary.csv',summary)]:a.write_csv(O/name,rows)
    result=dict(points=len(spatial),native_states=len(points),native_qualified=sum(r['qualified']=='True' for r in points),identity_passed=sum(r['current_identity_relative']<=1e-12 and r['field_identity_max_abs']==0 for r in spatial),observer_final_state_exact=all(r['observer_potential_max_V']==r['observer_n_max_relative']==r['observer_p_max_relative']==0 for r in spatial),observer_field_max_V_per_cm=max(r['observer_field_max_V_per_cm'] for r in spatial),cell_gradient_max_V_per_cm=max(r['cell_gradient_max_V_per_cm'] for r in spatial),prior_diagnostic_floor_max_relative=max(r['max_relative'] for r in summary if r['mode']=='prior_diagnostic_floor'),mathematical_e_max_relative=max(r['max_relative'] for r in summary if r['mode']=='mathematical_floor' and r['carrier']=='e'),mathematical_h_max_relative=max(r['max_relative'] for r in summary if r['mode']=='mathematical_floor' and r['carrier']=='h'),production_parameters_fitted=False,production_changed=False,acceptance_changed=False)
    a.write(O/'summary.json',result)
    d.matrix.freeze(O/'calibration_evidence.json',[Path(__file__).resolve(),O/'export_evidence.json',FLOORS,O/'summary.json']+[O/n for n in ('field_identity.csv','observer_samples.csv','cell_formation.csv','spatial_summary.csv','formation_summary.csv')])
    print(result,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('unpack','export','analyze'));p.add_argument('--sha');args=p.parse_args()
    if args.action=='unpack':unpack(args.sha)
    elif args.action=='export':export()
    else:analyze()
