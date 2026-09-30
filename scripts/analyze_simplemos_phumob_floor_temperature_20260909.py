"""Compare native G floors with exact, Charon and Genius algorithms at four T."""
import math
from functools import lru_cache
from pathlib import Path
from scipy.optimize import brentq
import simplemos_phumob_floor_temperature_20260909 as t
import analyze_simplemos_phumob_native_cells_20260909 as c

a,d=t.a,t.d;OUT=t.OUT/'analysis'


def G(P,m,T):
    return 1-.89233*(.41372+P*(T/300/m)**.28227)**(-.19778)+.005978*(P*(300/T*m)**.72169)**(-1.80618)


def dG(P,m,T):
    x=(T/300/m)**.28227;y=(300/T*m)**.72169
    return .89233*.19778*x*(.41372+P*x)**(-1.19778)-.005978*1.80618*y**(-1.80618)*P**(-2.80618)


@lru_cache(None)
def minimum(mode,m,T):
    if mode=='exact':return brentq(lambda P:dG(P,m,T),.001,10,xtol=1e-15)
    if mode=='genius':
        tableT=max(10,min(990,int(T/10)*10));lo,hi=.01,1.
        while hi-lo>1e-3:
            p1=lo+.33*(hi-lo);p2=hi-.33*(hi-lo)
            if G(p1,m,tableT)<G(p2,m,tableT):hi=p2
            else:lo=p1
        return (lo+hi)/2
    P=.3246;step=1.
    for iteration in range(500):
        if abs(step)<=1e-5:return P
        x=(T/300/m)**.28227;v=.005978*1.80618/(.89233*.19778)*(m/(T/300))**(-.72169*1.80618)
        f=v*P**(-2.80618)-x*(.41372+P*x)**(-1.19778)
        df=-2.80618*v*P**(-3.80618)+1.19778*x*x*(.41372+P*x)**(-2.19778)
        step=f/df;P-=step
    raise AssertionError('Charon iteration bound reached')


def main():
    a.verify(t.OUT/'export_evidence.json')
    c.s.LOCAL,c.s.OUT=t.LOCAL,t.OUT
    source_files=[Path('D:/code-repo/Genius-TCAD-Open/src/material/Si/Si_mob_Philips.cc'),Path('D:/code-repo/tcad-charon/src/evaluators/Charon_Mobility_PhilipsThomas_impl.hpp')]
    # External source identities are retained separately; the freeze helper
    # records repository-local inputs only.
    a.write(OUT/'source_identity.json',{str(p):a.sha(p) for p in source_files})
    a.write(OUT/'contract.json',dict(algorithms=['exact G derivative root','Charon Newton, correction 1e-5, 500 iterations','Genius 10 K table, ternary .33 interval, width 1e-3'],
        inferred_floor='One first fully clamped native cell per carrier and T; diagnostic only, never a production parameter.',cell_relative_gate=1e-7,temperature_gate_K=1e-8))
    d.matrix.freeze(OUT/'freeze.json',[Path(__file__).resolve(),OUT/'contract.json',OUT/'source_identity.json',t.OUT/'export_evidence.json'])
    summaries=[];floorrows=[];details=[];temperaturechecks=[]
    for job in a.rows(t.OUT/'native_points.csv'):
        T=float(job['temperature_K']);data=c.load(job);f=data['fields']
        exp=t.LOCAL/'exports/phumob'/job['key'];manifest=a.read(exp/'field_manifest.json')['fields']
        temp=next(r for r in manifest if 'Temperature' in r['name'] and r['region']==0)
        assert temp['unit']=='K'
        actual=c.scalar(exp/'fields'/f"{temp['name']}_region0.csv")
        error=max(abs(v-T) for v in actual.values());assert error<=1e-8
        temperaturechecks.append(dict(key=job['key'],prescribed_K=T,export_max_error_K=error))
        for car in ('e','h'):
            mx,mn,theta,nref,alpha,m,otherm=(1417,52.2,2.285,9.68e16,.68,1,1.258) if car=='e' else (470.5,44.9,2.247,2.23e17,.719,1.258,1)
            lattice=mx*(T/300)**(-theta);muN=mx*mx/(mx-mn)*(T/300)**(3*alpha-1.5);muC=mx*mn/(mx-mn)*(300/T)**.5
            nodes={}
            for node in f['n']:
                nd,na,n,p=[f[key][node] for key in ('nd','na','n','p')]
                nd*=1+nd*nd/(.21*nd*nd+4e20**2);na*=1+na*na/(.5*na*na+7.2e20**2)
                other=p if car=='e' else n;nsc=nd+na+other
                P=(T/300)**2/(2.459*nsc**(2/3)/3.97e13+3.828*(n+p)/(1.36e20*m))
                z=P**.6478;F=(.7643*z+2.2999+6.5502*m/otherm)/(z+2.3670-.8552*m/otherm)
                A=muN*nsc*(nref/nsc)**alpha+muC*(n+p)
                nodes[node]=(P,A,(nd if car=='e' else na)+other/F,na if car=='e' else nd)
            def predict(cid,mode,floor=None):
                vs=data['cells'][cid];measures=[float(v['measure_um2']) for v in vs]
                values=[]
                for v,w in zip(vs,measures):
                    P,A,C,I=nodes[data['mapping'][int(v['vertex'])]]
                    limit=minimum('exact' if mode=='inferred' else mode,m,T)
                    gg=floor if mode=='inferred' and P<limit else G(max(P,limit),m,T)
                    values.append(w*lattice*A/(A+lattice*(C+I*gg)))
                return math.fsum(values)/math.fsum(measures)
            train=next(cid for cid in sorted(data['cells']) if all(nodes[data['mapping'][int(v['vertex'])]][0]<minimum('exact',m,T) and float(v['measure_um2'])>0 for v in data['cells'][cid]))
            target=data['native'][car][train]
            inferred=brentq(lambda g:predict(train,'inferred',g)-target,0,1,xtol=1e-15)
            floorrows.append(dict(key=job['key'],temperature_K=T,carrier=car,training_cell=train,inferred_G_floor=inferred,
                exact_P=minimum('exact',m,T),exact_G=G(minimum('exact',m,T),m,T),
                genius_P=minimum('genius',m,T),genius_G=G(minimum('genius',m,T),m,T),
                charon_P=minimum('charon',m,T),charon_G=G(minimum('charon',m,T),m,T)))
            for mode in ('exact','charon','genius','inferred'):
                errors=[];unclamped=[]
                for cid in data['cells']:
                    relative=predict(cid,mode,inferred)/data['native'][car][cid]-1
                    if cid!=train:errors.append(abs(relative))
                    if all(nodes[data['mapping'][int(v['vertex'])]][0]>=minimum('exact',m,T) for v in data['cells'][cid]):unclamped.append(abs(relative))
                    details.append(dict(key=job['key'],carrier=car,mode=mode,cell=cid,relative=relative))
                summaries.append(dict(key=job['key'],carrier=car,mode=mode,heldout_max_relative=max(errors),unclamped_max_relative=max(unclamped),qualified=max(errors)<=1e-7))
    for name,rows in [('temperatures.csv',temperaturechecks),('floors.csv',floorrows),('cells.csv',details),('summary.csv',summaries)]:a.write_csv(OUT/name,rows)
    d.matrix.freeze(OUT/'evidence.json',[OUT/'freeze.json']+[OUT/n for n in ('temperatures.csv','floors.csv','cells.csv','summary.csv')])
    print(floorrows,flush=True)


if __name__=='__main__':main()
