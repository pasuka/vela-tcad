"""Held-out HFS replay with previously frozen PhuMob cutoff witnesses.

No fitting, no production parameter promotion. This isolates an old low-field
qualification gap from the new HFS transformation and averaging order.
"""
import argparse
import math
from pathlib import Path
import calibrate_simplemos_hfs_formation_v2_20260914 as f
import audit_simplemos_enormal_floor_control_20260912 as old

a,d,c=f.a,f.d,f.c
O=f.O/'prior_floor_control_20260914'


def main(stage):
    assert not (O/f'{stage}_evidence.json').exists()
    a.verify(f.OUT/f'{stage}_evidence.json');a.verify(c.O/'floor_position_evidence.json')
    roots={r['carrier']:r for r in a.rows(old.ROOTS) if r['key']=='T300'}
    data=[r for r in a.rows(f.OUT/f'{stage}_cells.csv') if r['mode']=='boundary_any_contact_all_cut1']
    results=[];summary=[]
    for job in f.native.jobs(stage):
        if job['arm']!='observed':continue
        geo=c.geometry(job['device']);fields=f.fields(f.L/f'{stage}_exports'/job['name'])
        en={}
        for cid,cell in geo['cells'].items():
            ns=cell['nodes'];psi=fields['psi']
            vec=f.np.array([psi[ns[1]]-psi[ns[0]],psi[ns[2]]-psi[ns[0]]])@geo['gradient'][cid]*1e4
            en[cid]=abs(float(vec@geo['gd'][cid]))
        for b,car in enumerate(('e','h')):
            for mode in ('P_lower','P_upper'):
                bulk={}
                for node in fields['n']:
                    state={q:fields[q][node] for q in ('nd','na','n','p')}
                    nd,na,n,p=(state[q] for q in ('nd','na','n','p'))
                    nd*=1+nd*nd/(.21*nd*nd+4e20**2);na*=1+na*na/(.5*na*na+7.2e20**2)
                    mass=(1,1.258)[b]
                    P=1/(2.459*(nd+na+(p if b==0 else n))**(2/3)/3.97e13+3.828*(n+p)/(1.36e20*mass))
                    G=float(roots[car]['G_floor']) if P<float(roots[car][mode]) else 1-.89233*(.41372+P*(1/mass)**.28227)**(-.19778)+.005978*(P*mass**.72169)**(-1.80618)
                    bulk[node]=c.f.mobility(state,car,G)
                rr=[]
                for row in data:
                    if row['key']!=job['key'] or row['carrier']!=car:continue
                    cid=int(row['cell']);nodes=geo['cells'][cid]['nodes'];F=float(row['F_V_cm'])
                    low=[1/(1/bulk[k]+c.inverse_surface(en[cid],geo['distance'][k],fields['nd'][k]+fields['na'][k],b)) for k in nodes]
                    weights=geo['weights'][cid]
                    value=math.fsum(w*f.canali(m,F,b) for w,m in zip(weights,low))
                    before=f.canali(math.fsum(w*m for w,m in zip(weights,low)),F,b)
                    er=value/float(row['native_mu'])-1
                    out=dict(key=job['key'],carrier=car,mode=mode,cell=cid,relative_error=er,
                        mean_before_saturation_relative=before/float(row['native_mu'])-1)
                    results.append(out);rr.append(out)
                summary.append(dict(key=job['key'],carrier=car,mode=mode,cells=len(rr),
                    max_relative=max(abs(r['relative_error']) for r in rr),
                    mean_first_max_relative=max(abs(r['mean_before_saturation_relative']) for r in rr)))
        print(stage,job['key'],'prior floor control',flush=True)
    a.write_csv(O/f'{stage}_cells.csv',results);a.write_csv(O/f'{stage}_summary.csv',summary)
    result=dict(new_fitted_parameters=0,production_parameter_qualified=False,
        left_max_relative=max(r['max_relative'] for r in summary if r['mode']=='P_lower'),
        right_max_relative=max(r['max_relative'] for r in summary if r['mode']=='P_upper'),
        scope='Held-out prior diagnostic floor; does not repair the original mathematical floor gate.')
    a.write(O/f'{stage}_summary.json',result)
    d.matrix.freeze(O/f'{stage}_evidence.json',[Path(__file__).resolve(),old.ROOTS,c.O/'floor_position_evidence.json',
        f.OUT/f'{stage}_evidence.json',O/f'{stage}_cells.csv',O/f'{stage}_summary.csv',O/f'{stage}_summary.json'])
    print(result,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--stage',choices=('pilot','rest'),required=True)
    main(p.parse_args().stage)
