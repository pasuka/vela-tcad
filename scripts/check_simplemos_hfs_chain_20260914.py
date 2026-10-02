"""Independent high-precision HFS element chain versus state finite differences.

Diagnostic prototype only: no C++ production Jacobian qualification. Includes
PhuMob's two population partials, all three Enormal potential columns, HFS
contact/boundary branches, vertex saturation and native box weights. The 1 V/cm
empirical branch is tested away from its boundary, never called differentiable
at the boundary. The original PhuMob mathematical floor remains unfitted.
"""
import argparse
import math
from decimal import Decimal as D, localcontext
from pathlib import Path
import numpy as np
import calibrate_simplemos_hfs_formation_v2_20260914 as f

a,d,c=f.a,f.d,f.c
O=f.O/'chain_20260914'
STEPS=('1e-22','5e-23')


def dec(x):return D(str(float(x)))
def dot(x,y):return sum((u*v for u,v in zip(x,y)),D(0))
def norm(x):return dot(x,x).sqrt()


def evaluate(state, geom, b):
    grad, projector, gd, weights, distances, impurities, baseline, contact = geom
    vt=dec(c.c.bgn.VT)
    gradient=[tuple(sum((state[k][v]*grad[v][j] for v in range(3)),D(0)) for j in range(2)) for k in range(3)]
    signed=dot(gradient[0],gd);en=abs(signed)
    direction=gradient[0] if contact else tuple(dot(row,gradient[1+b]) for row in projector)
    raw=norm(direction);F=raw if raw>=1 else D(0)
    # Fields and F are native units (V/cm). All state columns are volts.
    dF=[D(0)]*9
    if F:
        for j in range(3):
            v=grad[j] if contact else tuple(dot(row,grad[j]) for row in projector)
            dF[(0 if contact else 1+b)*3+j]=dot(direction,v)/raw
    den=[(D(1) if signed>0 else D(-1) if signed<0 else D(0))*dot(grad[j],gd) for j in range(3)]
    velocity,beta=(map(D,('1.07e7','1.109')) if b==0 else map(D,('8.37e6','1.213')))
    B,C,lam,delta,eta=[dec(x) for x in c.PARAMS[b]]
    result=D(0);jac=[D(0)]*9
    for j in range(3):
        ref=baseline[j]
        n=ref['n']*((state[0][j]-ref['psi']-state[1][j]+ref['phin'])/vt).exp()
        p=ref['p']*((state[2][j]-ref['phip']-state[0][j]+ref['psi'])/vt).exp()
        bulk,dp,_=c.h.hp((ref['nd'],ref['na'],n,p),b)
        damp=(-distances[j]/D('.01')).exp()
        cc=C*(impurities[j]+1)**lam
        inv=D(0);dinv=D(0)
        if en:
            z=cc*en**(D(2)/3);denom=B+z
            inv=damp*(en/denom+en*en/delta+en**3/eta)
            dinv=damp*((B+z/3)/(denom*denom)+2*en/delta+3*en*en/eta)
        low=1/(1/bulk+inv)
        dl=[-low*low*dinv*den[k] if k<3 else D(0) for k in range(9)]
        coeff=(low/bulk)**2/vt
        dl[j]+=coeff*(dp[0]-dp[1]);dl[3+j]-=coeff*dp[0];dl[6+j]+=coeff*dp[1]
        if F:
            z=low*F/velocity;A=1+z**beta
            mu=low*A**(-1/beta);dm=A**(-1-1/beta)
            df=-low*low/velocity*z**(beta-1)*dm
        else:mu=low;dm=D(1);df=D(0)
        result+=weights[j]*mu
        for k in range(9):jac[k]+=weights[j]*(dm*dl[k]+df*dF[k])
    return result,jac,raw,en


def main(stage):
    assert not (O/f'{stage}_evidence.json').exists()
    a.verify(f.OUT/f'{stage}_evidence.json')
    a.write(O/f'{stage}_contract.json',dict(precision=110,steps_V=STEPS,
        relative_gate=1e-12,zero_derivative_absolute_gate=1e-60,
        selection='Every cell adjoining 792/1000/1009/1057/1089; first cell in each contact-count and boundary-count class; closest raw field below/above 1 V/cm for each carrier.',
        state_columns='All 9 local psi,phin,phip columns; both carriers. Doping,BGN,temperature fixed; densities respond with the frozen thermal voltage.',
        scope='Independent analytic prototype versus Decimal state differences, not production matrix Jv; discontinuous cutoff location is not qualified.',
        production_changed=False,acceptance_changed=False))
    checks=[];selected=[]
    formation=a.rows(f.OUT/f'{stage}_cells.csv')
    with localcontext() as ctx:
        ctx.prec=110
        for job in f.native.jobs(stage):
            if job['arm']!='observed':continue
            g,contacts,boundary=f.geometry(job['device'])
            fields=f.fields(f.L/f'{stage}_exports'/job['name'])
            ids=set()
            for node in (792,1000,1009,1057,1089):ids.update(g['nodecells'].get(node,[]))
            classes=set()
            for cid,cell in sorted(g['cells'].items()):
                key=(sum(n in contacts for n in cell['nodes']),len(boundary[cid]))
                if key not in classes:ids.add(cid);classes.add(key)
            for car in ('e','h'):
                rr=[r for r in formation if r['key']==job['key'] and r['carrier']==car and r['mode']=='boundary_qf_any_contact']
                for side in (False,True):
                    xx=[r for r in rr if (float(r['F_V_cm'])>=1)==side]
                    if xx:ids.add(int(min(xx,key=lambda r:abs(float(r['F_V_cm'])-1))['cell']))
            for cid in sorted(ids):
                ns=g['cells'][cid]['nodes']
                st=g['gradient'][cid]*1e4
                gg=np.vstack((-st.sum(axis=0),st))
                grad=[tuple(dec(v) for v in row) for row in gg]
                # Ensure exact partition-of-unity after decimal conversion.
                grad[0]=tuple(-grad[1][i]-grad[2][i] for i in range(2))
                proj=np.eye(2)
                if boundary[cid]:
                    proj=boundary[cid][0][1]
                    if any(np.linalg.norm(p-proj)>1e-10 for _,p in boundary[cid][1:]):proj=np.zeros((2,2))
                refs=[{k:dec(fields[k][n]) for k in ('psi','phin','phip','n','p','nd','na')} for n in ns]
                state=[[r[k] for r in refs] for k in ('psi','phin','phip')]
                geom=(grad,[[dec(v) for v in row] for row in proj],tuple(dec(v) for v in g['gd'][cid]),
                    [dec(v) for v in g['weights'][cid]],[dec(g['distance'][n]) for n in ns],
                    [dec(fields['nd'][n]+fields['na'][n]) for n in ns],refs,any(n in contacts for n in ns))
                for b,car in enumerate(('e','h')):
                    value,jac,raw,en=evaluate(state,geom,b)
                    assert abs(raw-1)>D('1e-8'), 'Branch boundary requires a separate one-sided study'
                    selected.append(dict(key=job['key'],cell=cid,carrier=car,raw_F=float(raw),Enormal=float(en),
                        contact_vertices=sum(n in contacts for n in ns),boundary_edges=len(boundary[cid]),
                        cutoff_active=raw<1,mobility=float(value)))
                    for col in range(9):
                        slopes=[]
                        for sh in STEPS:
                            h=D(sh);xp=[x.copy() for x in state];xm=[x.copy() for x in state]
                            xp[col//3][col%3]+=h;xm[col//3][col%3]-=h
                            slopes.append((evaluate(xp,geom,b)[0]-evaluate(xm,geom,b)[0])/(2*h))
                        analytical=jac[col]
                        error=abs(slopes[-1]/analytical-1) if analytical else abs(slopes[-1])
                        step=abs((slopes[0]-slopes[1])/analytical) if analytical else abs(slopes[0]-slopes[1])
                        gate=D('1e-12') if analytical else D('1e-60')
                        checks.append(dict(key=job['key'],cell=cid,carrier=car,column=('psi','phin','phip')[col//3],node=ns[col%3],
                            analytic=float(analytical),fd=float(slopes[-1]),error=float(error),step_error=float(step),
                            exact_zero_analytic=not analytical,qualified=max(error,step)<=gate))
            print(stage,job['key'],'chain checks',len(checks),flush=True)
    a.write_csv(O/f'{stage}_checks.csv',checks);a.write_csv(O/f'{stage}_selection.csv',selected)
    summary=dict(element_carrier_pairs=len(selected),columns=len(checks),qualified=sum(r['qualified'] for r in checks),
        max_nonzero_relative=max((r['error'] for r in checks if not r['exact_zero_analytic']),default=0),
        max_zero_absolute=max((r['error'] for r in checks if r['exact_zero_analytic']),default=0),
        all_prototype_columns_qualified=all(r['qualified'] for r in checks),
        native_matrix_qualified=False,production_matrix_qualified=False,cutoff_boundary_qualified=False)
    a.write(O/f'{stage}_summary.json',summary)
    d.matrix.freeze(O/f'{stage}_evidence.json',[Path(__file__).resolve(),Path(c.h.__file__),f.OUT/f'{stage}_evidence.json',
        O/f'{stage}_contract.json',O/f'{stage}_checks.csv',O/f'{stage}_selection.csv',O/f'{stage}_summary.json'])
    print(summary,flush=True);assert summary['all_prototype_columns_qualified']


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--stage',choices=('pilot','rest'),required=True)
    main(p.parse_args().stage)
