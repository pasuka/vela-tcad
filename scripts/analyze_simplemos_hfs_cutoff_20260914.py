"""Qualify Load/Plot integrity before interpreting the native HFS cutoff probe."""
import math
from pathlib import Path
from collections import defaultdict
from scipy.spatial import cKDTree
import probe_simplemos_hfs_native_cutoff_v3_20260914 as p

a,d,c,L,O=p.a,p.d,p.c,p.L,p.O


def main():
    assert not (O/'analysis_evidence.json').exists()
    a.verify(O/'export_evidence.json')
    requested=defaultdict(dict)
    for r in a.rows(O/'requested.csv'):requested[int(r['index'])][int(r['node'])]=r
    integrity=[];identities=[];loaded={}
    for i,F in enumerate(p.FORCES):
        for arm in ('baseline','observed'):
            src=L/'exports'/arm/f's{i:02d}'
            fields={key:c.scalar(src,name+'_region0.csv') for key,name in
                [('psi','ElectrostaticPotential'),('phin','eQuasiFermiPotential'),('phip','hQuasiFermiPotential'),('n','eDensity'),('p','hDensity')]}
            errors={key:max(abs(values[k]-float(requested[i][k][key])) if key in ('psi','phin','phip') else abs(values[k]/float(requested[i][k][key])-1) for k in values) for key,values in fields.items()}
            good=errors['psi']<=1e-12 and max(errors['phin'],errors['phip'])<=1e-10 and max(errors['n'],errors['p'])<=1e-8
            integrity.append(dict(index=i,F_requested=F,arm=arm,**errors,qualified=good));loaded[i,arm]=fields
        rows=c.identity(L/'exports/baseline'/f's{i:02d}',L/'exports/observed'/f's{i:02d}')
        identities+= [dict(index=i,**r) for r in rows]
    a.write_csv(O/'load_integrity.csv',integrity);a.write_csv(O/'observer_identity.csv',identities)
    integrity_summary=dict(all_loaded_states_qualified=all(r['qualified'] for r in integrity),
        max_identity_relative=max(r['max_relative'] for r in identities),
        max_identity_absolute=max(r['max_absolute'] for r in identities),
        observer_identity_qualified=all(r['max_relative']<=1e-12 for r in identities))
    a.write(O/'integrity_summary.json',integrity_summary)
    d.matrix.freeze(O/'integrity_evidence.json',[Path(__file__).resolve(),O/'export_evidence.json',O/'requested.csv',
        O/'load_integrity.csv',O/'observer_identity.csv',O/'integrity_summary.json'])
    print(integrity_summary,flush=True)
    assert integrity_summary['all_loaded_states_qualified'] and integrity_summary['observer_identity_qualified']
    xy={int(r['id']):(float(r['x_um']),float(r['y_um'])) for r in a.rows(L/'exports/observed/s00/nodes.csv')}
    ids=sorted(xy);tree=cKDTree([xy[k] for k in ids])
    rows=[];excluded=[]
    for car in ('e','h'):
        for r in a.rows(L/'raw/bundle/observed'/f'hfs_observer_{car}_Silicon_1.csv'):
            distance,j=tree.query([float(r['x_um']),float(r['y_um'])]);node=ids[int(j)];assert distance<1e-12
            candidates=[]
            for i in range(len(p.FORCES)):
                z=loaded[i,'observed']
                mismatch=max(abs(float(r['n_cm3'])/z['n'][node]-1),abs(float(r['p_cm3'])/z['p'][node]-1))
                if mismatch<=1e-13:candidates.append((i,mismatch))
            if len(candidates)!=1:
                excluded.append(dict(carrier=car,node=node,sequence=r['last_call'],matches=len(candidates)))
                continue
            i,mismatch=candidates[0];F=float(r['F_V_cm']);m=float(r['mulow_cm2_V_s']);mu=float(r['mu_cm2_V_s'])
            expected=0. if p.FORCES[i]<1 else p.FORCES[i]
            # The exactly-one test is a boundary witness, not a smooth branch gate.
            good=(F==0 if p.FORCES[i]<1 else abs(F-p.FORCES[i])<=1e-7) if p.FORCES[i]!=1 else (F==0 or abs(F-1)<=1e-7)
            rows.append(dict(index=i,carrier=car,node=node,sequence=r['last_call'],F_requested=p.FORCES[i],
                F_received=F,state_mismatch=mismatch,mulow=m,mu=mu,saturation_ratio=mu/m-1,
                field_error=abs(F-expected),threshold_boundary=p.FORCES[i]==1,qualified=good))
    summaries=[]
    for i,F in enumerate(p.FORCES):
        for car in ('e','h'):
            rr=[r for r in rows if r['index']==i and r['carrier']==car]
            assert rr,'No unambiguously identified native samples for state'
            summaries.append(dict(index=i,carrier=car,F_requested=F,samples=len(rr),nodes=len({r['node'] for r in rr}),
                min_received=min(r['F_received'] for r in rr),max_received=max(r['F_received'] for r in rr),
                zero_received=sum(r['F_received']==0 for r in rr),qualified=sum(r['qualified'] for r in rr),
                saturation_ratio_min=min(r['saturation_ratio'] for r in rr),saturation_ratio_max=max(r['saturation_ratio'] for r in rr)))
    a.write_csv(O/'samples.csv',rows);a.write_csv(O/'excluded_calls.csv',excluded);a.write_csv(O/'threshold_summary.csv',summaries)
    result=dict(samples=len(rows),excluded_calls=len(excluded),all_samples_qualified=all(r['qualified'] for r in rows),
        native_threshold_bracket_V_cm=[.999999,1.000001],
        exact_boundary_differentiable=False,DC_qualified=False,
        scope='This fixed-state native probe isolates field formation; contact-QF and boundary projection disabled only here. Does not qualify a self-consistent response or production Jacobian.')
    a.write(O/'summary.json',result)
    d.matrix.freeze(O/'analysis_evidence.json',[Path(__file__).resolve(),O/'integrity_evidence.json',
        O/'samples.csv',O/'excluded_calls.csv',O/'threshold_summary.csv',O/'summary.json'])
    print(result,summaries,flush=True);assert result['all_samples_qualified']


if __name__=='__main__':main()
