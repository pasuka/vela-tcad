"""Held-out Enormal states against previously frozen left/right cutoff controls."""
import math
import calibrate_simplemos_enormal_20260912 as p
a,d=p.a,p.d
ROOTS=p.R/'reference_tcad/simplemos_sentaurus2022/phumob_floor_20260909/ramped/cutoff/roots.csv'

def run():
    a.verify(p.O/'calibration_evidence.json')
    roots={r['carrier']:r for r in a.rows(ROOTS) if r['key']=='T300'}
    values={}
    for row in a.rows(p.O/'cell_formation.csv'):
        if row['mode']=='mathematical_floor':values[row['key'],row['carrier'],int(row['cell'])]=row
    rows=[];summary=[]
    for job in a.read(p.O/'native_contract.json')['jobs']:
        if job['arm']!='observed':continue
        key=job['name'].removesuffix('_observed');src=p.L/'exports'/job['name'];geo=p.geometry(job['device'])
        fields={k:p.scalar(src,v+'_region0.csv') for k,v in [('nd','DonorConcentration'),('na','AcceptorConcentration'),('n','eDensity'),('p','hDensity')]}
        for b,car in enumerate(('e','h')):
            for mode in ('P_lower','P_upper'):
                mu={};interval=0
                for k in fields['n']:
                    state={q:fields[q][k] for q in fields};nd,na,n,h=(state[q] for q in ('nd','na','n','p'))
                    nd*=1+nd*nd/(.21*nd*nd+4e20**2);na*=1+na*na/(.5*na*na+7.2e20**2)
                    mass=(1,1.258)[b];P=1/(2.459*(nd+na+(h if b==0 else n))**(2/3)/3.97e13+3.828*(n+h)/(1.36e20*mass))
                    lo=float(roots[car]['P_lower']);hi=float(roots[car]['P_upper']);interval+=lo<P<hi
                    G=float(roots[car]['G_floor']) if P<float(roots[car][mode]) else 1-.89233*(.41372+P*(1/mass)**.28227)**(-.19778)+.005978*(P*mass**.72169)**(-1.80618)
                    mu[k]=p.f.mobility(state,car,G)
                errs=[]
                for cid,cell in geo['cells'].items():
                    ref=values[key,car,cid];F=float(ref['field_V_per_cm']);target=float(ref['native_mu'])
                    vals=[1/(1/mu[k]+p.inverse_surface(F,geo['distance'][k],fields['nd'][k]+fields['na'][k],b)) for k in cell['nodes']]
                    trial=math.fsum(w*v for w,v in zip(geo['weights'][cid],vals));err=trial/target-1;errs.append(abs(err))
                    rows.append(dict(key=key,carrier=car,mode=mode,cell=cid,relative_error=err))
                summary.append(dict(key=key,carrier=car,mode=mode,interval_nodes=interval,max_relative=max(errs),below_original_cell_gate=max(errs)<=1e-7))
    a.write_csv(p.O/'floor_position_control.csv',summary);a.write_csv(p.O/'floor_position_cells.csv',rows)
    result=dict(new_fitted_parameters=0,production_changed=False,left_max_relative=max(r['max_relative'] for r in summary if r['mode']=='P_lower'),right_max_relative=max(r['max_relative'] for r in summary if r['mode']=='P_upper'))
    a.write(p.O/'floor_position_summary.json',result)
    d.matrix.freeze(p.O/'floor_position_evidence.json',[p.O/'calibration_evidence.json',ROOTS,p.Path(__file__).resolve()]+[p.O/n for n in ('floor_position_control.csv','floor_position_cells.csv','floor_position_summary.json')])
    print(result,flush=True)
if __name__=='__main__':run()
