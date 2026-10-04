"""Independent Decimal constitutive/flux reference for the unified state contract."""
from decimal import Decimal as Q,localcontext
from pathlib import Path
import math
import validate_simplemos_split_contract_v2_20260911 as v
F=lambda x:Q.from_float(float(x))

def reference(data,state,precision):
    with localcontext() as ctx:
        ctx.prec=precision;N=len(data['nodes']);a=[[F(z) for z in n] for n in data['nodes']]
        scale=F(state['potential_scale_V']);x=[F(hi)+F(lo) for hi,lo in state['coordinates']]
        psi=[x[i]*scale for i in range(N)];en=[x[N+i]*scale+F(state['electron_reference_V'][i]) for i in range(N)];hp=[x[2*N+i]*scale+F(state['hole_reference_V'][i]) for i in range(N)]
        c=[F(z) for z in data['phumob']['common']];params=[[F(z) for z in data['phumob'][k]] for k in ('electron','hole')]
        def raw(p,b):return 1-c[10]*(c[11]+p*(1/c[4+b])**c[13])**(-c[15])+c[12]*(p*c[4+b]**c[14])**(-c[16])
        roots=[]
        for b in range(2):
            lo,hi=Q(-80),Q(80);first=(1/c[4+b])**c[13];second=c[4+b]**c[14]
            for _ in range(260):
                mid=(lo+hi)/2;p=mid.exp();pos=c[10]*c[15]*first*(c[11]+p*first)**(-c[15]-1);neg=c[12]*c[16]*(p*second)**(-c[16])/p
                if pos>neg:hi=mid
                else:lo=mid
            p=((lo+hi)/2).exp();roots.append((p,abs(raw(p,b))))
        def mu(b,nd,na,n,p):
            def cluster(z,ref,g):u=z/ref;return z*(1+u*u/(1+g*u*u))
            nd=cluster(nd*c[17],c[0]*c[17],c[2]);na=cluster(na*c[17],c[1]*c[17],c[3]);n*=c[17];p*=c[17]
            m=params[b];mx=m[0]*c[18];mn=m[1]*c[18];sc=nd+na+(n if b else p)
            if not sc:return mx/c[18]
            mass=c[4+b];ratio=mass/c[5-b];free=n+p;screen=1/(c[6]*sc**(Q(2)/3)/F(3.97e13)+c[7]*free/(F(1.36e20)*mass))
            power=screen**F(.6478);f=(F(.7643)*power+F(2.2999)+F(6.5502)*ratio)/(power+F(2.3670)-F(.8552)*ratio)
            g=roots[b][1] if screen<roots[b][0] else max(roots[b][1],raw(screen,b));eff=na+g*nd+c[9]*n/f if b else nd+g*na+c[8]*p/f
            mus=mx*mx/(mx-mn)*(sc/eff)*(m[3]*c[17]/sc)**m[4]+mx*mn/(mx-mn)*free/eff
            return 1/(1/mx+1/mus)/c[18]
        n=[Q(0)]*N;p=[Q(0)]*N;mn=[Q(0)]*N;mp=[Q(0)]*N;srh=[Q(0)]*N;r=[Q(0)]*(3*N);magnitude=[Q(0)]*(3*N)
        for i,aa in enumerate(a):
            if aa[9]<=0:continue
            n[i]=aa[9]*((psi[i]-en[i])/aa[10]).exp();p[i]=aa[9]*((hp[i]-psi[i])/aa[10]).exp();phys=[F(z) for z in data['node_physics'][i]]
            mn[i]=mu(0,phys[0],phys[1],n[i],p[i]);mp[i]=mu(1,phys[0],phys[1],n[i],p[i])
            # Deliberately use the unfactored density product for the independent SRH check.
            srh[i]=Q(0) if en[i]==hp[i] else (n[i]*p[i]-aa[9]**2)/(phys[3]*(n[i]+aa[9])+phys[2]*(p[i]+aa[9]));source=srh[i]*phys[4]
            r[N+i]+=source;r[2*N+i]+=source;magnitude[N+i]+=abs(source);magnitude[2*N+i]+=abs(source)
        for i,j,g in data['poisson_edges']:
            f=F(g)*(psi[i]-psi[j]);r[i]+=f;r[j]-=f;magnitude[i]+=abs(f);magnitude[j]+=abs(f)
        for i,aa in enumerate(a):
            qn=aa[15]*aa[16]*n[i]*aa[12];qp=aa[15]*aa[16]*p[i]*aa[13];qd=aa[15]*aa[16]*aa[11]*aa[14]
            r[i]=(r[i]+qn-qp-qd-aa[18])/aa[17];magnitude[i]=(magnitude[i]+abs(qn)+abs(qp)+abs(qd)+abs(aa[18]))/abs(aa[17])
        def B(z):return Q(1) if z==0 else z/(z.exp()-1)
        edges=[];ports={};conducting=set()
        for e in data['transport_edges']:
            i,j=e['i'],e['j'];fe=fh=Q(0)
            if e['weights']:
                me=sum(F(w)*mn[k] for k,w in e['weights']);mh=sum(F(w)*mp[k] for k,w in e['weights']);vt=a[i][10]
                eta=(psi[j]-psi[i])/vt+(a[j][9]/a[i][9]).ln();etaH=(psi[j]-psi[i])/vt+(a[i][9]/a[j][9]).ln()
                # Density SG, independent from the C++ factored expm1 form.
                fe=Q(0) if en[i]==en[j] else me*F(e['coefficient'])*(n[i]*B(-eta)-n[j]*B(eta));fh=Q(0) if hp[i]==hp[j] else mh*F(e['coefficient'])*(p[i]*B(etaH)-p[j]*B(-etaH));conducting.update((i,j))
            edges.append((fe,fh));r[N+i]+=fe;r[N+j]-=fe;r[2*N+i]+=fh;r[2*N+j]-=fh
            magnitude[N+i]+=abs(fe);magnitude[N+j]+=abs(fe);magnitude[2*N+i]+=abs(fh);magnitude[2*N+j]+=abs(fh)
            for name,sign in e['ports'].items():ports[name]=ports.get(name,Q(0))+F(sign)*(fe-fh)*F(data['current_factor'])
        for i in range(N):
            if i not in conducting:r[N+i]=en[i]/scale;r[2*N+i]=hp[i]/scale;magnitude[N+i]=magnitude[2*N+i]=Q(1)
        for b,bc in enumerate(data['boundary']):
            for i,z in bc:r[b*N+i]=x[b*N+i]-F(z);magnitude[b*N+i]=Q(1)
        return dict(nodes=[(n[i],p[i],mn[i],mp[i],srh[i],r[i],r[N+i],r[2*N+i],psi[i]) for i in range(N)],edges=edges,ports=ports,scale=magnitude)

def main():
    v.verify(v.OUT/'run_evidence.json');records=[];directions=[];checks=[]
    for c in v.rows(v.OUT/'runs.csv'):
        dest=Path(c['dest']);data=v.read(dest/'operator.json');N=len(data['nodes']);rr=v.rows(dest/'results.csv');group={}
        for r in rr:group.setdefault(r['case'],{})[r['kind'],r['index']]=tuple(r[k] for k in 'abcdefghi')
        assert group['base']==group['partition']
        base=v.rows(dest/'baseline.csv');model_error=0.;residual_error=Q(0);worst_row=-1;worst_old=Q(0);worst_new=Q(0)
        expected_base=None;cpp_base=None
        for label in ('base','micro_0','micro_1','micro_2'):
            state=v.read(dest/('results.csv.'+label+'.json'));h=reference(data,state,160);l=reference(data,state,120)
            max_ref=Q(0);agreement=Q(0);values=[];pred=[]
            for i in range(N):
                cpp=[Q(z) for z in group[label]['node',str(i)]];values.extend(cpp);pred.extend(h['nodes'][i])
                for k,(actual,truth) in enumerate(zip(cpp,h['nodes'][i])):
                    scale=max(abs(truth),Q('1e-300')) if k not in (5,6,7) else max(h['scale'][(k-5)*N+i],Q('1e-300'))
                    max_ref=max(max_ref,abs(actual-truth)/scale);agreement=max(agreement,abs(l['nodes'][i][k]-truth)/scale)
                if label=='base':
                    for k,col in ((0,'n'),(1,'p')):
                        prod=F(data['baseline_density'][i][k])
                        if prod:model_error=max(model_error,float(abs(cpp[k]/prod-1)))
                    for k,col in enumerate(('psi_residual','phin_residual','phip_residual')):
                        if k==0 or data['nodes'][i][9]>0:
                            err=abs(cpp[k+5]-F(base[i][col]))/max(h['scale'][k*N+i],Q('1e-300'))
                            if err>residual_error:residual_error=err;worst_row=k*N+i;worst_old=F(base[i][col]);worst_new=cpp[k+5]
            for ei,e in enumerate(h['edges']):
                cpp=group[label]['edge',str(ei)]
                for k in range(2):
                    scale=max(abs(e[k]),Q('1e-200'));max_ref=max(max_ref,abs(Q(cpp[k])-e[k])/scale)
            for name,p in h['ports'].items():max_ref=max(max_ref,abs(Q(group[label]['port',name][0])-p)/max(abs(p),Q('1e-200')))
            assert max_ref<Q('1e-60') and agreement<Q('1e-70'),(c,label,max_ref,agreement)
            if label=='base':expected_base=pred;cpp_base=values
            else:
                with localcontext() as ctx:
                    ctx.prec=160
                    diff=[x-y for x,y in zip(values,cpp_base)];ref=[x-y for x,y in zip(pred,expected_base)]
                    # Compare each changed observable relatively; ignore mathematically zero responses.
                    active=[(x,y) for x,y in zip(diff,ref) if y!=0 and abs(y)>Q('1e-85')*max(abs(x),Q(1))]
                    response=max((abs(x-y)/abs(y) for x,y in active),default=Q(0));assert active and response<Q('1e-10'),(c,label,response)
                checks.append(dict(device=c['device'],vd=c['vd'],index=c['index'],label=label,changed=sum(x!=0 for x in diff),max_response_relative=float(response),reference_relative=float(max_ref)))
            print(c['device'],c['vd'],c['index'],label,'reference',float(max_ref),flush=True)
        # Compare live element/edge mobilities against the original production values.
        for e in data['transport_edges']:
            for k in range(2):
                value=sum(F(w)*Q(group['base']['node',str(i)][2+k]) for i,w in e['weights'])
                old=F(e['baseline_mu'][k])
                if old:model_error=max(model_error,float(abs(value/old-1)))
        assert model_error<2e-12,(c,model_error)
        # Retain the frozen reconstruction failure while completing independent checks.
        print('production residual',c['device'],c['vd'],c['index'],float(residual_error),worst_row,flush=True)
        dr=v.rows(Path(c['old_dest'])/'trial_0.csv')
        for size,alpha in (('full',Q(1)),('small',Q(1)/4096)):
            for b in range(3):
                with localcontext() as ctx:
                    ctx.prec=120
                    fd=[(Q(group['plus_'+size]['node',str(i)][5+b])-Q(group['minus_'+size]['node',str(i)][5+b]))/(2*alpha) for i in range(N)]
                    jv=[F(dr[b*N+i]['JdxR'])-F(dr[b*N+i]['R']) for i in range(N)]
                    err=sum((x-y)**2 for x,y in zip(fd,jv)).sqrt()/max(sum(x*x for x in fd).sqrt(),sum(x*x for x in jv).sqrt(),Q('1e-300'))
                directions.append(dict(device=c['device'],vd=c['vd'],index=c['index'],block=b,size=size,relative=float(err),qualified=err<Q('1e-4')))
        records.append(dict(device=c['device'],vd=c['vd'],index=c['index'],partition=True,checkpoint_jobs=9,max_model_relative=model_error,max_scaled_residual_difference=float(residual_error),production_residual_qualified=residual_error<Q('1e-10'),worst_row=worst_row,worst_old=str(worst_old),worst_new=str(worst_new)))
    for name,data in (('summary.csv',records),('micro.csv',checks),('directions.csv',directions)):v.csvout(v.OUT/'analysis_v3'/name,data)
    v.freeze(v.OUT/'analysis_v3/analysis_evidence.json',[Path(__file__).resolve(),v.OUT/'run_evidence.json']+[v.OUT/'analysis_v3'/n for n in ('summary.csv','micro.csv','directions.csv')])
    print(records,flush=True)

if __name__=='__main__':main()
