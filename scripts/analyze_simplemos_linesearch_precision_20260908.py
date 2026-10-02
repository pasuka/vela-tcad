"""Independent Decimal residual/merit scan of actual and ideal Newton trials."""
import argparse
from pathlib import Path
from decimal import Decimal as D, localcontext, getcontext
import math
import audit_simplemos_linesearch_precision_20260908 as audit
import analyze_simplemos_minority_residual_20260906 as hp

a=audit.a; d=audit.d; OUT=audit.OUT/'precision'; MODES=('rounded','kernel','split','ideal')


def prepare():
    a.verify(audit.OUT/'freeze.json')
    a.write(OUT/'contract.json',dict(modes=dict(
        rounded='100-digit sum of exported rounded production edge fluxes and source products; isolates summation/norm rounding.',
        kernel='Exact binary G/psi/n/p Poisson operands and exported SG/SRH kernel operands, fixed double material/geometry coefficients.',
        split='Exact x times potential scale, retained qf reference plus increment, recomputed densities and SG/SRH in high precision at actual rounded candidate x.',
        ideal='As split, but x=base_x+alpha*capped_step formed without double rounding; diagnostic mathematical trial, not a representable production iterate.'),
        precisions=[60,100],gates=dict(precision_agreement_over_row_scale=1e-20,rounded_assembly_reconstruction_over_scale=1e-10),
        decisions='Use each representation own baseline merit. Reproduce actual norm and global-mode formula; norm<base or norm<=(1-1e-4*alpha)*base. Never accept a state from a diagnostic scan.',
        scope='All Poisson rows and both free Si carrier rows; exact binary operands, all source and boundary/gauge rows retained. Linear defect uses exported full Jacobian and raw/capped directions.'))
    d.matrix.freeze(OUT/'freeze.json',[Path(__file__).resolve(),Path(hp.__file__),audit.OUT/'freeze.json',OUT/'contract.json'])


def records(path):
    return [{k:hp.exact(x) for k,x in r.items()} for r in a.rows(path)]


def merit(R,pn,meta):
    N=len(pn);block=[sum((x*x for x in R[b*N:(b+1)*N]),D(0)) for b in range(3)]
    if meta['residual_mode']=='l2':value=sum(block)
    else:value=sum(block[b]*hp.exact(meta['weight_'+f])/hp.exact(meta['scale_'+f])**2 for b,f in enumerate(('psi','n','p')))
    if meta['global_mode']!='off':
        for b,f in ((1,'n'),(2,'p')):
            source=sum((R[b*N+i] for i in range(N) if not pn[i]['bc'+f]),D(0))
            value+=(source/hp.exact(meta['global_e_scale' if b==1 else 'global_h_scale']))**2
    return value.sqrt(),[b.sqrt() for b in block]


def evaluate(data,base,precision):
    with localcontext() as ctx:
        ctx.prec=precision;pn,edges,nodes,pe,dr,meta=data;N=len(pn);result={}
        prod=[r['R'] for r in dr]
        result['production']=dict(R=prod,scale=[max(abs(x),D('1e-300')) for x in prod],merit=merit(prod,pn,meta)[0])
        for mode in MODES:
            ns=[dict(n) for n in nodes];ps=[dict(n) for n in pn]
            if mode in ('split','ideal'):
                for i in range(N):
                    for b,f in enumerate(('xpsi','xn','xp')):
                        xx=ps[i][f] if mode=='split' else base[b*N+i]['x']+hp.exact(meta['alpha'])*base[b*N+i]['capped_step']
                        ps[i][f]=xx
                    ns[i]['psi']=ps[i]['xpsi']*ps[i]['potential_scale']
                    ns[i]['einc']=ps[i]['xn']*ps[i]['potential_scale'];ns[i]['hinc']=ps[i]['xp']*ps[i]['potential_scale']
                    ps[i]['psi']=ns[i]['psi'];ni=ns[i]['ni'];vt=ns[i]['Vt']
                    if ni>0:
                        ps[i]['n']=ni*hp.clamp((ns[i]['psi']-ns[i]['eref']-ns[i]['einc'])/vt).exp()
                        ps[i]['p']=ni*hp.clamp((ns[i]['href']+ns[i]['hinc']-ns[i]['psi'])/vt).exp()
                    else:ps[i]['n']=ps[i]['p']=D(0)
            R=[D(0) for _ in range(3*N)];scale=[D(0) for _ in R]
            for e in pe:
                i,j=int(e['node0']),int(e['node1']);flux=e['rounded_flux'] if mode=='rounded' else e['G']*(ps[i]['psi']-ps[j]['psi'])
                R[i]+=flux;R[j]-=flux;scale[i]+=abs(flux);scale[j]+=abs(flux)
            for i,n in enumerate(ps):
                if mode=='rounded':sn=n['rounded_electron_charge'];sp=n['rounded_hole_charge'];sd=n['rounded_dopant_charge']
                else:
                    factor=n['q']*n['area_factor'];sn=factor*n['n']*n['voln'];sp=factor*n['p']*n['volp'];sd=factor*n['net_doping']*n['vold']
                R[i]=(R[i]+sn-sp-sd-n['interface_rhs'])/n['poisson_scale']
                scale[i]=(scale[i]+abs(sn)+abs(sp)+abs(sd)+abs(n['interface_rhs']))/n['poisson_scale']
            for b,car,field in ((1,'electron','nflux'),(2,'hole','pflux')):
                for e in edges:
                    i,j=int(e['node0']),int(e['node1'])
                    f=e[field] if mode=='rounded' else hp.edge_flux(e,ns,car,'kernel' if mode=='kernel' else 'partition')
                    R[b*N+i]+=f;R[b*N+j]-=f;scale[b*N+i]+=abs(f);scale[b*N+j]+=abs(f)
                for i,n in enumerate(ns):
                    source=n['rate']*n['volume']*n['source_factor']/n['scale'] if mode=='rounded' else hp.srh(n,'kernel' if mode=='kernel' else 'partition')
                    R[b*N+i]+=source;scale[b*N+i]=max(scale[b*N+i],abs(source))
            for i in range(N):
                for b,f in enumerate(('psi','n','p')):
                    k=b*N+i
                    if pn[i]['bc'+f] or (b and nodes[i]['ni']<=0):
                        R[k]=prod[k]
                        if mode=='ideal':R[k]+=ps[i][('xpsi','xn','xp')[b]]-pn[i][('xpsi','xn','xp')[b]]
                        scale[k]=max(abs(R[k]),D(1))
                    scale[k]=max(scale[k],D('1e-300'))
            m,blocks=merit(R,pn,meta);result[mode]=dict(R=R,scale=scale,merit=m,blocks=blocks)
        return result


def run():
    getcontext().prec=120;a.verify(OUT/'freeze.json');a.verify(audit.OUT/'replay_evidence.json')
    rows=[];worst=[];linears=[];recon=[];coordinates=[];files=[OUT/'freeze.json',audit.OUT/'replay_evidence.json']
    for job in a.read(audit.OUT/'contract.json')['jobs']:
        root=Path(job['dest']);base=records(root/'base_direction.csv');node=job['node'];N=len(base)//3
        J=records(root/'jacobian.csv')
        for label in ('raw_step','capped_step'):
            res=[r['R'] for r in base]
            for e in J:res[int(e['row'])]+=e['value']*base[int(e['column'])][label]
            linears.append(dict(key=job['case']['key'],direction=label,raw_equals_capped=all(r['raw_step']==r['capped_step'] for r in base),
                defect_l2=float(sum(x*x for x in res).sqrt()),worst_hole_defect=float(res[2*N+node]),worst_hole_linear_ratio=float(abs(res[2*N+node])/max(abs(base[2*N+node]['R']),D('1e-300')))))
        baseline=None
        for label in ['base']+['trial_'+str(k) for k in range(13)]:
            prefix=root/label;pn=records(str(prefix)+'_poisson_nodes.csv');ed=records(str(prefix)+'_edges.csv');nd=records(str(prefix)+'_nodes.csv');pe=records(str(prefix)+'_poisson_edges.csv');dr=records(str(prefix)+'_direction.csv');meta=a.rows(str(prefix)+'_merit.csv')[0]
            data=(pn,ed,nd,pe,dr,meta);evaluations={prec:evaluate(data,base,prec) for prec in (60,100)};high=evaluations[100]
            if label=='base':baseline=high
            agreement=max(abs(evaluations[60][m]['R'][i]-high[m]['R'][i])/high[m]['scale'][i] for m in MODES for i in range(3*N))
            reconstruction=max(abs(high['rounded']['R'][i]-dr[i]['R'])/high['rounded']['scale'][i] for i in range(3*N))
            assert agreement<=D('1e-20'),(label,agreement)
            assert reconstruction<=D('1e-10'),(label,reconstruction)
            normerr=abs(float(high['production']['merit'])/float(meta['merit'])-1)
            assert normerr<=1e-12,(label,normerr)
            recon.append(dict(key=job['case']['key'],label=label,precision_agreement=float(agreement),rounded_reconstruction=float(reconstruction),norm_replay_relative=normerr))
            alpha=hp.exact(meta['alpha'])
            for mode in ('production',)+MODES:
                e=high[mode];current=baseline[mode]['merit'];m=e['merit'];rr=e['R'][2*N+node]
                passes=label!='base' and (m<current or m<=(D(1)-D('1e-4')*alpha)*current)
                rows.append(dict(key=job['case']['key'],label=label,alpha=float(alpha),mode=mode,merit=str(m),relative_change=float(m/current-1),would_decrease=passes,
                    hole_node=node,hole_residual=str(rr),hole_row_ratio=float(abs(rr)/e['scale'][2*N+node]) if mode!='production' else '',
                    psi_norm=float(e['blocks'][0]) if mode!='production' else '',n_norm=float(e['blocks'][1]) if mode!='production' else '',p_norm=float(e['blocks'][2]) if mode!='production' else ''))
            for b,f in enumerate(('psi','phin','phip')):
                expected=[alpha*base[b*N+i]['capped_step'] for i in range(N)]
                actual=[dr[b*N+i]['state']-base[b*N+i]['x'] for i in range(N)]
                coordinates.append(dict(key=job['case']['key'],label=label,field=f,nonzero_intended=sum(x!=0 for x in expected),nonzero_actual=sum(x!=0 for x in actual),
                    max_rounding_error_V=float(max(abs(x-y) for x,y in zip(actual,expected))*pn[0]['potential_scale'])))
            if label in ('base','trial_0'):
                order=sorted(range(N),key=lambda i:abs(high['kernel']['R'][i]-dr[i]['R']),reverse=True)[:10]
                for i in order:worst.append(dict(key=job['case']['key'],label=label,node=i,production_Rpsi=float(dr[i]['R']),rounded_Rpsi=float(high['rounded']['R'][i]),kernel_Rpsi=float(high['kernel']['R'][i]),split_Rpsi=float(high['split']['R'][i]),ideal_Rpsi=float(high['ideal']['R'][i])))
            print(job['case']['key'],label,'recon',float(reconstruction),'production',rows[-5]['relative_change'],'split',rows[-2]['relative_change'],flush=True)
    for name,data in [('merit',rows),('checks',recon),('linear_defect',linears),('coordinates',coordinates),('worst_poisson',worst)]:a.write_csv(OUT/(name+'.csv'),data)
    d.matrix.freeze(OUT/'evidence.json',files+list(OUT.glob('*.csv')))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('prepare','run'));globals()[parser.parse_args().action]()
