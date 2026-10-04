"""Independent split-coordinate Poisson response and mixed-path merit gates."""
import decimal
from pathlib import Path
import validate_simplemos_split_psi_v2_20260911 as v
D=decimal.Decimal.from_float
Z=decimal.Decimal(0)

def evaluate(data,trial,kernel,mode,precision):
    with decimal.localcontext() as ctx:
        ctx.prec=precision
        nodes=[[D(float(x)) for x in a] for a in data['nodes']];n=len(nodes)
        psi=[(D(float(r['hi']))+D(float(r['lo'])))*a[6] for r,a in zip(kernel,nodes)]
        result=[Z]*n;scale=[Z]*n
        for ei,ej,g in data['edges']:
            f=D(float(g))*(psi[ei]-psi[ej]);result[ei]+=f;result[ej]-=f;scale[ei]+=abs(f);scale[ej]+=abs(f)
        for i,a in enumerate(nodes):
            ne=ho=Z
            if a[9]>0:
                xn=D(float(trial['candidate'][n+i]));xp=D(float(trial['candidate'][2*n+i]))
                ne=a[9]*max(decimal.Decimal(-500),min(decimal.Decimal(500),(psi[i]-a[7]-xn*a[6])/a[10])).exp()
                ho=a[9]*max(decimal.Decimal(-500),min(decimal.Decimal(500),(a[8]+xp*a[6]-psi[i])/a[10])).exp()
            qn=a[15]*a[16]*ne*a[12];qp=a[15]*a[16]*ho*a[13];qd=a[15]*a[16]*a[11]*a[14]
            result[i]=(result[i]+qn-qp-qd-a[18])/a[17]
            scale[i]=(scale[i]+abs(qn)+abs(qp)+abs(qd)+abs(a[18]))/abs(a[17])
        for node,value in data['bcpsi']:
            result[node]=(D(float(kernel[node]['hi']))+D(float(kernel[node]['lo'])))-D(float(value));scale[node]=decimal.Decimal(1)
        return result,scale

def main():
    v.verify(v.OUT/'replay_evidence.json');details=[];checks=[];summaries=[]
    for c in v.rows(v.OUT/'identity.csv'):
        dest=Path(c['new_dest']);old=Path(c['dest']);data=v.read(dest/'inputs.json');n=len(data['nodes'])
        raw=v.rows(dest/'kernel.csv');group={}
        for r in raw:group.setdefault((r['trial'],r['mode']),[]).append(r)
        starts={};own=[];max_kernel=Z;max_agreement=Z;exact_coordinates=0;rounded_id=0
        for trial in data['trials']:
            label=trial['label'];dr=v.rows(old/('trial_0.csv' if label=='base' else label+'.csv'))
            for mode in ('rounded','split'):
                kr=group[label,mode]
                low,sc=evaluate(data,trial,kr,mode,60);rr,ss=evaluate(data,trial,kr,mode,100)
                agreement=max(abs(x-y)/max(s,decimal.Decimal('1e-300')) for x,y,s in zip(low,rr,ss))
                err=max(abs(decimal.Decimal(r['R'])-x)/max(s,decimal.Decimal('1e-300')) for r,x,s in zip(kr,rr,ss))
                assert agreement<decimal.Decimal('1e-25') and err<decimal.Decimal('1e-25'),(c,label,mode,agreement,err)
                max_kernel=max(max_kernel,err);max_agreement=max(max_agreement,agreement)
                if mode=='split':
                    with decimal.localcontext() as ctx:
                        ctx.prec=120
                        for i,r in enumerate(kr):
                            assert D(float(r['hi']))+D(float(r['lo']))==D(trial['x'][i])+D(trial['alpha'])*D(trial['step'][i])
                            exact_coordinates+=1
                if mode=='rounded':
                    field='R' if label=='base' else 'Rtrial'
                    assert all(float(r['R'])==float(dr[i][field]) for i,r in enumerate(kr)),(c,label,'baseline changed')
                    rounded_id+=n
                if label=='base':starts[mode]=rr
                with decimal.localcontext() as ctx:
                    ctx.prec=100
                    delta=sum((x-y)*(x+y)*D(float(dr[i]['weight']))/D(float(dr[i]['scale']))**2 for i,(x,y) in enumerate(zip(rr,starts[mode])))
                    carrier_delta=sum((D(float(r['Rtrial']))-D(float(r['R'])))*(D(float(r['Rtrial']))+D(float(r['R'])))*D(float(r['weight']))/D(float(r['scale']))**2 for r in dr[n:]) if label!='base' else Z
                    # For this isolated candidate only Poisson consumes the tail.
                    # Keep actual carrier residuals, so this is not full split DD.
                    coupled_delta=delta+carrier_delta
                    a=D(trial['alpha']);jv=[D(float(r['JdxR']))-D(float(r['R'])) for r in dr[:n]]
                    linear_error=sum((x-y-a*j)**2 for x,y,j in zip(rr,starts[mode],jv)).sqrt()
                    linear_scale=max(sum((a*j)**2 for j in jv).sqrt(),decimal.Decimal('1e-300'))
                rec=dict(device=c['device'],vd=c['vd'],index=c['index'],trial=label,alpha=trial['alpha'],mode=mode,poisson_square_delta=str(delta),held_carrier_square_delta=str(carrier_delta),mixed_square_delta=str(coupled_delta),mixed_decreases=label!='base' and coupled_delta<0,poisson_linear_relative=float(linear_error/linear_scale) if label!='base' else 0.)
                details.append(rec);own.append(rec)
        summary=dict(device=c['device'],vd=c['vd'],index=c['index'],rounded_decreases=sum(r['mixed_decreases'] for r in own if r['mode']=='rounded'),split_decreases=sum(r['mixed_decreases'] for r in own if r['mode']=='split'),trials=13,exact_split_coordinates=exact_coordinates,identical_rounded_residuals=rounded_id,max_kernel_scaled_error=str(max_kernel),max_decimal_agreement=str(max_agreement))
        summaries.append(summary);print(summary,flush=True)
    v.csvout(v.OUT/'comparison.csv',details);v.csvout(v.OUT/'summary.csv',summaries)
    v.freeze(v.OUT/'analysis_evidence.json',[Path(__file__).resolve(),v.OUT/'replay_evidence.json',v.OUT/'comparison.csv',v.OUT/'summary.csv'])

if __name__=='__main__':main()
