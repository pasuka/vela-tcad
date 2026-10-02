"""Independent BGN SG/SRH and Poisson reconstruction of actual backtracks."""
from decimal import Decimal as D, localcontext, getcontext
from pathlib import Path
import audit_simplemos_e1089_linesearch_20260909 as v
import analyze_simplemos_linesearch_precision_20260908 as core
import audit_simplemos_bgn_psi_precision_20260908 as bgn

s=v.s;hp=core.hp;OUT=v.OUT/'precision'

def edge_flux(e,nodes,carrier,mode):
    electron=carrier=='electron';prefix='e' if electron else 'h';coef=e['ncoef' if electron else 'pcoef']
    if coef==0:return D(0)
    i,j=int(e['node0']),int(e['node1']);vt=e['Vt'];ni0,ni1=e['ni0'],e['ni1'];assert ni0>0 and ni1>0
    if mode=='kernel':psi0,psi1,q0,q1=(e[prefix+k] for k in ('psi0','psi1','phi0','phi1'))
    else:
        p,q=nodes[i],nodes[j];psi0,psi1=p['psi'],q['psi'];q0=p[prefix+'ref']+p[prefix+'inc'];q1=q[prefix+'ref']+q[prefix+'inc']
    eta=(psi1-psi0)/vt+(ni1/ni0 if electron else ni0/ni1).ln()
    z0=(psi0-q0)/vt if electron else (q0-psi0)/vt;z1=(psi1-q1)/vt if electron else (q1-psi1)/vt
    b0,b1=hp.clamp(z0),hp.clamp(z1)
    drive=((q1-q0) if electron else (q0-q1))/vt+(b0-z0)-(b1-z1)
    return coef*hp.bern(eta if electron else -eta)*ni1*b1.exp()*hp.em1(drive)/e['scale']

def records(path):return [{k:hp.exact(x) for k,x in r.items()} for r in s.a.rows(Path(path))]

def run():
    s.a.verify(v.OUT/'replay_evidence.json');assert all(r['qualified']=='True' for r in s.a.rows(v.OUT/'identity.csv'))
    s.a.write(OUT/'contract.json',dict(scope='Actual first-iteration 13 rejected candidates from the frozen e1089 failure; all Poisson and free Si carrier rows. Output-only replay preserves state, status, repeated residuals.',precisions=[60,100],modes=core.MODES,BGN='Include intrinsic-density log drift in both carrier fluxes; independent unfactorized density-SG check on kernel operands.',coefficients='Production double geometry, ni, mobility and lifetimes are fixed per candidate. Split/ideal arithmetic is a fixed-coefficient diagnostic, not full high-precision constitutive or self-consistent recomputation.',gates=dict(precision_agreement=1e-20,reconstruction=1e-10,row=1e-6),decisions='Compare each representation against its own base. Preserve actual rejected candidates; no high precision diagnostic becomes an accepted state.'))
    s.d.matrix.freeze(OUT/'freeze.json',[Path(__file__).resolve(),Path(core.__file__),Path(hp.__file__),Path(bgn.__file__),v.OUT/'replay_evidence.json',OUT/'contract.json'])
    hp.edge_flux=edge_flux;getcontext().prec=120;root=v.LOCAL/'failed_e1089/trace';base=records(root/'base_direction.csv');N=len(base)//3;node=1089
    rows=[];checks=[];coordinates=[];hotspots=[];linear=[];coeff_checks=[];baseline=None
    J=records(root/'jacobian.csv')
    for name in ('raw_step','capped_step'):
        defect=[r['R'] for r in base]
        for e in J:defect[int(e['row'])]+=e['value']*base[int(e['column'])][name]
        linear.append(dict(direction=name,raw_equals_capped=all(r['raw_step']==r['capped_step'] for r in base),linear_e1089= str(defect[N+node]),linear_e1089_over_residual=float(abs(defect[N+node]/base[N+node]['R'])),linear_defect_norm=str(sum(x*x for x in defect).sqrt())))
    for label in ['base']+['trial_'+str(k) for k in range(13)]:
        pn,ed,nd,pe,dr=[records(root/(label+suffix+'.csv')) for suffix in ('_poisson_nodes','_edges','_nodes','_poisson_edges','_direction')];meta=s.a.rows(root/(label+'_merit.csv'))[0]
        data=(pn,ed,nd,pe,dr,meta);evaluations={precision:core.evaluate(data,base,precision) for precision in (60,100)};high=evaluations[100]
        if label=='base':baseline=high
        agreement=max(abs(evaluations[60][m]['R'][i]-high[m]['R'][i])/high[m]['scale'][i] for m in core.MODES for i in range(3*N))
        reconstruction=max(abs(high['rounded']['R'][i]-dr[i]['R'])/high['rounded']['scale'][i] for i in range(3*N))
        assert agreement<=D('1e-20'),(label,agreement);assert reconstruction<=D('1e-10'),(label,reconstruction)
        assert abs(float(high['production']['merit'])/float(meta['merit'])-1)<=1e-12
        # Check the factored reference against separately expressed density SG.
        with localcontext() as ctx:
            ctx.prec=100
            for car,prefix in (('electron','e'),('hole','h')):
                selected=[e for e in ed if e['ncoef' if car=='electron' else 'pcoef']!=0 and node in (int(e['node0']),int(e['node1']))]
                for e in selected:
                    values=[e[k] for k in ('ni0','ni1',prefix+'psi0',prefix+'psi1',prefix+'phi0',prefix+'phi1','Vt')]
                    expected=bgn.density_flux(values,car=='electron')*e['ncoef' if car=='electron' else 'pcoef']/e['scale']
                    observed=edge_flux(e,nd,car,'kernel');relative=abs(observed-expected)/max(abs(observed),abs(expected),D('1e-300'));assert relative<D('1e-40'),(label,car,relative)
                    coeff_checks.append(dict(label=label,carrier=car,edge=int(e['edge']),density_vs_factored_relative=float(relative)))
        checks.append(dict(label=label,precision_agreement=float(agreement),reconstruction=float(reconstruction)))
        alpha=hp.exact(meta['alpha'])
        for mode in ('production',)+core.MODES:
            e=high[mode];norm=e['merit'];start=baseline[mode]['merit'];rr=e['R'][N+node]
            delta=sum((y-x)*(y+x) for x,y in zip(baseline[mode]['R'],e['R']))
            scale=e['scale'][N+node] if mode!='production' else high['rounded']['scale'][N+node]
            rows.append(dict(label=label,alpha=float(alpha),mode=mode,merit=str(norm),norm_squared_delta=str(delta),relative_change=float(norm/start-1),would_decrease=label!='base' and delta<0,electron_residual=str(rr),electron_row_ratio=float(abs(rr)/scale),psi_norm=float(e['blocks'][0]) if mode!='production' else '',electron_norm=float(e['blocks'][1]) if mode!='production' else '',hole_norm=float(e['blocks'][2]) if mode!='production' else ''))
        for b,f in enumerate(('psi','phin','phip')):
            expected=[alpha*base[b*N+i]['capped_step'] for i in range(N)];actual=[dr[b*N+i]['state']-base[b*N+i]['x'] for i in range(N)]
            coordinates.append(dict(label=label,field=f,nonzero_intended=sum(x!=0 for x in expected),nonzero_actual=sum(x!=0 for x in actual),maximum_intended_V=float(max(abs(x) for x in expected)*pn[0]['potential_scale']),maximum_rounding_error_V=float(max(abs(x-y) for x,y in zip(actual,expected))*pn[0]['potential_scale'])))
        for mode in core.MODES:
            for i in sorted(range(N),key=lambda i:abs(high[mode]['R'][i]-dr[i]['R']),reverse=True)[:8]:hotspots.append(dict(label=label,mode=mode,node=i,production_Rpsi=str(dr[i]['R']),reference_Rpsi=str(high[mode]['R'][i]),difference=str(high[mode]['R'][i]-dr[i]['R'])))
        print(label,'e1089',rows[-5]['electron_row_ratio'],'production delta',rows[-5]['norm_squared_delta'],'ideal decrease',rows[-1]['would_decrease'],flush=True)
    for name,values in (('merit',rows),('checks',checks),('coordinates',coordinates),('poisson_hotspots',hotspots),('linear',linear),('independent_flux_checks',coeff_checks)):s.a.write_csv(OUT/(name+'.csv'),values)
    s.d.matrix.freeze(OUT/'evidence.json',[OUT/'freeze.json']+list(OUT.glob('*.csv')))

if __name__=='__main__':run()
