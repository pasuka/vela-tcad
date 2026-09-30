"""Independent source-aware residual reconstruction at actual rejected trials."""
from decimal import Decimal as D, getcontext
from pathlib import Path
import audit_simplemos_lowvd_source_precision_20260910 as v
import analyze_simplemos_e1089_precision_20260909 as previous

s=v.s;core=previous.core;hp=core.hp;OUT=v.OUT/'precision'

def run():
    s.a.verify(v.OUT/'replay_evidence.json');getcontext().prec=120
    s.a.write(OUT/'contract.json',dict(precisions=[60,100],source='Retain the exact binary alpha and eight frozen normalized pair-source values in both carrier equations. Rounded mode uses the actual combined SRH/source rate; kernel/split/ideal add the independently specified constant source.',scope='All rows at base and trial 0, 5, 12 of each of four final rejected searches; exact production squared-norm differences and update rounding for all 13 trials.',gates=dict(precision=1e-20,reconstruction=1e-10),coefficients='Fixed production double coefficients at each actual candidate; split and ideal are diagnostics, not accepted self-consistent states.'))
    s.d.matrix.freeze(OUT/'freeze.json',[Path(__file__).resolve(),Path(previous.__file__),Path(core.__file__),Path(hp.__file__),v.OUT/'replay_evidence.json',OUT/'contract.json'])
    hp.edge_flux=previous.edge_flux;original_srh=hp.srh
    merits=[];checks=[];coordinates=[];hotspots=[];linears=[];actual=[]
    for job in s.a.read(v.OUT/'contract.json')['jobs']:
        root=Path(job['dest']);base=previous.records(root/'base_direction.csv');N=len(base)//3
        lines=Path(job['source_file']).read_text().splitlines();assert int(lines[0].split()[0])==N
        source={int(i):hp.exact(value)*hp.exact(job['alpha']) for i,value in (line.split() for line in lines[1:])};assert len(source)==8
        hp.srh=lambda node,mode:original_srh(node,mode)+source.get(int(node['node']),D(0))
        J=previous.records(root/'jacobian.csv');defect=[r['R'] for r in base]
        for e in J:defect[int(e['row'])]+=e['value']*base[int(e['column'])]['raw_step']
        linears.append(dict(key=job['key'],label=job['label'],raw_equals_capped=all(r['raw_step']==r['capped_step'] for r in base),linear_defect_norm=str(sum(x*x for x in defect).sqrt()),base_norm=str(sum(r['R']**2 for r in base).sqrt())))
        baseline=None
        for label in ['base']+['trial_'+str(k) for k in range(13)]:
            dr=previous.records(root/(label+'_direction.csv'));meta=s.a.rows(root/(label+'_merit.csv'))[0];alpha=hp.exact(meta['alpha'])
            blocks=[sum((dr[b*N+i]['R']-base[b*N+i]['R'])*(dr[b*N+i]['R']+base[b*N+i]['R']) for i in range(N)) for b in range(3)]
            actual.append(dict(key=job['key'],source_label=job['label'],trial=label,alpha=float(alpha),delta_psi=str(blocks[0]),delta_n=str(blocks[1]),delta_p=str(blocks[2]),delta_total=str(sum(blocks)),decreases=sum(blocks)<0))
            for b,field in enumerate(('psi','phin','phip')):
                expected=[alpha*base[b*N+i]['capped_step'] for i in range(N)];observed=[dr[b*N+i]['state']-base[b*N+i]['x'] for i in range(N)]
                coordinates.append(dict(key=job['key'],source_label=job['label'],trial=label,field=field,intended_nonzero=sum(x!=0 for x in expected),actual_nonzero=sum(x!=0 for x in observed),max_scaled_rounding=float(max(abs(x-y) for x,y in zip(expected,observed)))))
            if label not in ('base','trial_0','trial_5','trial_12'):continue
            pn,ed,nd,pe=[previous.records(root/(label+suffix+'.csv')) for suffix in ('_poisson_nodes','_edges','_nodes','_poisson_edges')]
            data=(pn,ed,nd,pe,dr,meta);low=core.evaluate(data,base,60);high=core.evaluate(data,base,100)
            if label=='base':baseline=high
            agreement=max(abs(low[m]['R'][i]-high[m]['R'][i])/high[m]['scale'][i] for m in core.MODES for i in range(3*N))
            reconstruction=max(abs(high['rounded']['R'][i]-dr[i]['R'])/high['rounded']['scale'][i] for i in range(3*N))
            assert agreement<D('1e-20') and reconstruction<D('1e-10'),(job['key'],label,agreement,reconstruction)
            checks.append(dict(key=job['key'],source_label=job['label'],trial=label,precision_agreement=float(agreement),reconstruction=float(reconstruction)))
            for mode in ('production',)+core.MODES:
                e=high[mode];delta=sum((y-x)*(y+x) for x,y in zip(baseline[mode]['R'],e['R']))
                rowscale=high['rounded']['scale'] if mode=='production' else e['scale']
                free=[i for b,f in ((1,'n'),(2,'p')) for i in range(b*N,(b+1)*N) if not pn[i%N]['bc'+f] and nd[i%N]['ni']>0]
                worst=max(free,key=lambda i:abs(e['R'][i])/rowscale[i])
                merits.append(dict(key=job['key'],source_label=job['label'],trial=label,mode=mode,delta=str(delta),decreases=delta<0,merit=str(e['merit']),max_carrier_ratio=float(abs(e['R'][worst])/rowscale[worst]),worst_node=worst%N,worst_carrier='n' if worst<2*N else 'p',psi_norm=str(sum(x*x for x in e['R'][:N]).sqrt())))
                if mode!='production':
                    for i in sorted(range(N),key=lambda i:abs(e['R'][i]-dr[i]['R']),reverse=True)[:6]:
                        hotspots.append(dict(key=job['key'],source_label=job['label'],trial=label,mode=mode,node=i,production_Rpsi=str(dr[i]['R']),reference_Rpsi=str(e['R'][i])))
            print(job['key'],job['label'],label,'ratios',[(r['mode'],r['max_carrier_ratio'],r['decreases']) for r in merits[-5:]],flush=True)
    hp.srh=original_srh
    for name,data in (('merit',merits),('checks',checks),('coordinates',coordinates),('poisson_hotspots',hotspots),('linear',linears),('actual_deltas',actual)):s.a.write_csv(OUT/(name+'.csv'),data)
    s.d.matrix.freeze(OUT/'evidence.json',[OUT/'freeze.json']+list(OUT.glob('*.csv')))

if __name__=='__main__':run()
