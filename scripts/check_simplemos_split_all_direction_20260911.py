"""Separate QF increment rounding from the remaining Poisson directional defect."""
from decimal import Decimal as Q,localcontext
from pathlib import Path
import validate_simplemos_split_psi_v2_20260911 as v
F=lambda x:Q.from_float(float(x))

def residual(data,t,ideal,precision):
    with localcontext() as ctx:
        ctx.prec=precision;n=len(data['nodes']);a=[[F(z) for z in x] for x in data['nodes']]
        x=[F(b)+F(t['alpha'])*F(s) if ideal or i<n else F(t['candidate'][i]) for i,(b,s) in enumerate(zip(t['x'],t['step']))]
        psi=[x[i]*a[i][6] for i in range(n)];r=[Q(0)]*n
        for i,j,g in data['edges']:
            f=F(g)*(psi[i]-psi[j]);r[i]+=f;r[j]-=f
        for i,b in enumerate(a):
            ne=ho=Q(0)
            if b[9]>0:
                ne=b[9]*max(Q(-500),min(Q(500),(psi[i]-b[7]-x[n+i]*b[6])/b[10])).exp()
                ho=b[9]*max(Q(-500),min(Q(500),(b[8]+x[2*n+i]*b[6]-psi[i])/b[10])).exp()
            r[i]=(r[i]+b[15]*b[16]*(ne*b[12]-ho*b[13]-b[11]*b[14])-b[18])/b[17]
        for i,z in data['bcpsi']:r[i]=x[i]-F(z)
        return r

def main():
    v.verify(v.OUT/'analysis_evidence.json');records=[]
    for c in v.rows(v.OUT/'identity.csv'):
        data=v.read(Path(c['new_dest'])/'inputs.json');dr=v.rows(Path(c['dest'])/'trial_0.csv');n=len(data['nodes']);base=residual(data,data['trials'][0],True,100)
        for t in (data['trials'][1],data['trials'][-1]):
            for ideal in (False,True):
                with localcontext() as ctx:
                    ctx.prec=100
                    rr=residual(data,t,ideal,100);low=residual(data,t,ideal,60)
                    jv=[F(t['alpha'])*(F(r['JdxR'])-F(r['R'])) for r in dr[:n]]
                    norm=sum(z*z for z in jv).sqrt();defect=sum((x-y-z)**2 for x,y,z in zip(rr,base,jv)).sqrt()/norm
                    agree=sum((x-y)**2 for x,y in zip(low,rr)).sqrt()/norm;assert agree<Q('1e-20')
                    qf_error=max(abs(F(t['candidate'][i])-(F(t['x'][i])+F(t['alpha'])*F(t['step'][i]))) for i in range(n,3*n))
                    rec=dict(device=c['device'],vd=c['vd'],index=c['index'],alpha=t['alpha'],mode='all_coordinates' if ideal else 'psi_only',poisson_jv_relative=float(defect),precision_agreement=float(agree),max_normalized_qf_roundoff=float(qf_error),passes_previous_direction_gate=defect<Q('1e-4'))
                    records.append(rec);print(rec,flush=True)
    v.csvout(v.OUT/'all_direction.csv',records);v.freeze(v.OUT/'all_direction_evidence.json',[Path(__file__).resolve(),v.OUT/'analysis_evidence.json',v.OUT/'all_direction.csv'])

if __name__=='__main__':main()
