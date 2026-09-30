"""Distinguish clipping P from clipping G near the native PhuMob cutoff."""
import math
from pathlib import Path
from scipy.optimize import brentq
import simplemos_phumob_floor_ramp_20260909 as r
import analyze_simplemos_phumob_native_cells_20260909 as c
from analyze_simplemos_phumob_floor_temperature_20260909 import G,minimum

a,d=r.a,r.d;OUT=r.OUT/'cutoff'


def main():
    a.verify(r.OUT/'analysis/evidence.json')
    a.write(OUT/'contract.json',dict(scope='Use the already inferred floor to calculate both roots of G(P)=G_floor. Test P clipping at each root and G-value clipping on all other cells.',
        interpretation='Diagnostic identification only; root ambiguity is retained when the mesh does not sample the interval. No production parameters fitted.',
        original_cell_gate=1e-7,algorithm_distinction_relative=1e-12))
    d.matrix.freeze(OUT/'freeze.json',[Path(__file__).resolve(),r.OUT/'analysis/evidence.json',OUT/'contract.json'])
    floors=a.rows(r.OUT/'analysis/floors.csv');c.s.LOCAL,c.s.OUT=r.LOCAL,r.OUT
    results=[];cellrows=[];roots=[]
    for job in a.rows(r.OUT/'native_points.csv'):
        T=float(job['temperature_K']);data=c.load(job);f=data['fields']
        for car in ('e','h'):
            fit=next(x for x in floors if x['key']==job['key'] and x['carrier']==car)
            floor=float(fit['inferred_G_floor']);train=int(fit['training_cell'])
            mx,mn,theta,nref,alpha,m,otherm=(1417,52.2,2.285,9.68e16,.68,1,1.258) if car=='e' else (470.5,44.9,2.247,2.23e17,.719,1.258,1)
            exact=minimum('exact',m,T)
            low=brentq(lambda P:G(P,m,T)-floor,exact*.5,exact,xtol=1e-15)
            high=brentq(lambda P:G(P,m,T)-floor,exact,exact*2,xtol=1e-15)
            lattice=mx*(T/300)**(-theta);muN=mx*mx/(mx-mn)*(T/300)**(3*alpha-1.5);muC=mx*mn/(mx-mn)*(300/T)**.5
            values={mode:{} for mode in ('P_lower','P_upper','G_floor')};witnesses=0
            for node in f['n']:
                nd,na,n,p=[f[k][node] for k in ('nd','na','n','p')]
                nd*=1+nd*nd/(.21*nd*nd+4e20**2);na*=1+na*na/(.5*na*na+7.2e20**2)
                other=p if car=='e' else n;nsc=nd+na+other
                P=(T/300)**2/(2.459*nsc**(2/3)/3.97e13+3.828*(n+p)/(1.36e20*m))
                witnesses+=low<P<high
                z=P**.6478;F=(.7643*z+2.2999+6.5502*m/otherm)/(z+2.3670-.8552*m/otherm)
                A=muN*nsc*(nref/nsc)**alpha+muC*(n+p)
                for mode,gg in [('P_lower',G(max(P,low),m,T)),('P_upper',G(max(P,high),m,T)),('G_floor',floor if P<exact else max(G(P,m,T),floor))]:
                    eff=(nd+gg*na if car=='e' else na+gg*nd)+other/F
                    values[mode][node]=lattice*A/(A+lattice*eff)
            roots.append(dict(key=job['key'],carrier=car,temperature_K=T,G_floor=floor,P_exact=exact,P_lower=low,P_upper=high,interval_vertices=witnesses))
            for mode,nodal in values.items():
                errors=[]
                for cid,vs in data['cells'].items():
                    weights=[float(v['measure_um2']) for v in vs]
                    mu=math.fsum(w*nodal[data['mapping'][int(v['vertex'])]] for v,w in zip(vs,weights))/math.fsum(weights)
                    error=mu/data['native'][car][cid]-1
                    if cid!=train:errors.append(abs(error))
                    cellrows.append(dict(key=job['key'],carrier=car,cell=cid,mode=mode,relative=error))
                results.append(dict(key=job['key'],carrier=car,mode=mode,heldout_max_relative=max(errors),cell_gate_passed=max(errors)<=1e-7,distinction_gate_passed=max(errors)<=1e-12))
    a.write_csv(OUT/'roots.csv',roots);a.write_csv(OUT/'results.csv',results);a.write_csv(OUT/'cells.csv',cellrows)
    d.matrix.freeze(OUT/'evidence.json',[OUT/'freeze.json',OUT/'roots.csv',OUT/'results.csv',OUT/'cells.csv'])
    print(roots,results,flush=True)


if __name__=='__main__':main()
