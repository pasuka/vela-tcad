"""Isolated PhuMob carrier derivatives: Decimal chain rule versus independent FD.

No production model changes. Positive concentrations in cm^-3; mobility in
cm^2/(V s). The 300 K equations are from the local T-2022.03 manual, 270-283.
"""
from decimal import Decimal as D,localcontext
from functools import lru_cache
from pathlib import Path
import argparse,math,subprocess
import simplemos_phumob_calibration_native_20260908 as p

a,d=p.a,p.d
LOCAL=p.LOCAL/'cross_derivatives';OUT=p.OUT/'cross_derivatives'
NODES=(792,1000,1009,1057)
STEPS=(1e-4,3e-5,1e-5,3e-6,1e-6,3e-7)


@lru_cache(None)
def gminimum(carrier):
    m=D('1') if carrier==0 else D('1.258')
    t=(1/m)**D('.28227');u=m**D('.72169')
    def dg(x):return D('.89233')*D('.19778')*t*(D('.41372')+x*t)**D('-1.19778')-D('.005978')*D('1.80618')*u**D('-1.80618')*x**D('-2.80618')
    low,high=D('1e-12'),D('100')
    assert dg(low)<0 and dg(high)>0
    for _ in range(240):
        mid=(low+high)/2
        if dg(mid)>0:high=mid
        else:low=mid
    x=(low+high)/2
    value=1-D('.89233')*(D('.41372')+x*t)**D('-.19778')+D('.005978')*(x*u)**D('-1.80618')
    return x,abs(value)


def hp(state,carrier):
    nd,na,n,h=state
    nd=nd*(1+nd*nd/(D('.21')*nd*nd+D('4e20')**2))
    na=na*(1+na*na/(D('.5')*na*na+D('7.2e20')**2))
    mu_max,mu_min,nref,alpha,m,other_m=(map(D,('1417','52.2','9.68e16','.68','1','1.258')) if carrier==0 else map(D,('470.5','44.9','2.23e17','.719','1.258','1')))
    other=h if carrier==0 else n
    dnsc=(D(0),h) if carrier==0 else (n,D(0))
    free=n+h;nsc=nd+na+other
    aa=D('2.459')/D('3.97e13');bb=D('3.828')/(D('1.36e20')*m)
    screening=1/(aa*nsc**(D(2)/3)+bb*free)
    dp=[-screening**2*(aa*(D(2)/3)*nsc**(-D(1)/3)*x+bb*y) for x,y in zip(dnsc,(n,h))]
    z=screening**D('.6478');num=D('.7643')*z+D('2.2999')+D('6.5502')*m/other_m
    den=z+D('2.3670')-D('.8552')*m/other_m
    ff=num/den
    dff=D('.6478')*screening**D('-.3522')*(D('.7643')*den-num)/(den*den)
    t=(1/m)**D('.28227');u=m**D('.72169')
    floor_x,floor_g=gminimum(carrier)
    if screening<floor_x:
        gg,dgg=floor_g,D(0)
    else:
        raw=1-D('.89233')*(D('.41372')+screening*t)**D('-.19778')+D('.005978')*(screening*u)**D('-1.80618')
        gg=max(raw,floor_g)
        dgg=(D('.89233')*D('.19778')*t*(D('.41372')+screening*t)**D('-1.19778')-D('.005978')*D('1.80618')*u**D('-1.80618')*screening**D('-2.80618')) if raw>floor_g else D(0)
    minority_impurity=na if carrier==0 else nd
    neff=(nd if carrier==0 else na)+gg*minority_impurity+other/ff
    dneff=[(dgg*minority_impurity-other*dff/(ff*ff))*x+y/ff for x,y in zip(dp,dnsc)]
    mu_n=mu_max*mu_max/(mu_max-mu_min);mu_c=mu_max*mu_min/(mu_max-mu_min)
    first=mu_n*nref**alpha*nsc**(1-alpha)
    scattering=(first+mu_c*free)/neff
    dsc=[(first*(1-alpha)/nsc*x+mu_c*y-scattering*z)/neff for x,y,z in zip(dnsc,(n,h),dneff)]
    mu=mu_max*scattering/(mu_max+scattering)
    deriv=[mu_max**2/(mu_max+scattering)**2*x for x in dsc]
    return mu,deriv,screening<floor_x


def prepare():
    q=p.qualified
    a.verify(q.OUT/'final_evidence.json');a.verify(q.OUT/'comparison_evidence.json')
    samples=[];files=[Path(__file__).resolve(),p.LOCAL/'phumob_scalar_probe.exe',p.REPO/'scripts/diagnostics/simplemos_phumob_scalar_probe.cpp',p.REPO/'src/physics/MobilityModel.cpp',p.REPO/'include/vela/physics/MobilityModel.h',q.OUT/'final_evidence.json']
    for c in a.rows(q.OUT/'comparison.csv'):
        if c['model']!='old_slotboom' or int(c['index']) not in p.INDICES:continue
        state=q.LOCAL/'dc/old_slotboom'/c['case']/'vela'/f"vg_{int(c['index']):03d}"/'attempt_0/state.csv'
        values={int(x['node_id']):x for x in a.rows(state)}
        fields=p.b.LOCAL/'native_exports/old_slotboom'/c['case']/f"vg_{int(c['index']):03d}"/'fields'
        nd={int(x['node_id']):float(x['component0']) for x in a.rows(fields/'DonorConcentration_region0.csv')}
        na={int(x['node_id']):float(x['component0']) for x in a.rows(fields/'AcceptorConcentration_region0.csv')}
        for node in NODES:
            r=values[node]
            samples.append(dict(id=str(len(samples)),case=c['case'],index=int(c['index']),node=node,nd=nd[node],na=na[node],n=float(r['electrons_m3'])*1e-6,p=float(r['holes_m3'])*1e-6))
        files += [state,fields/'DonorConcentration_region0.csv',fields/'AcceptorConcentration_region0.csv']
    a.write(OUT/'contract.json',dict(samples=samples,precision=100,reference_steps_log_density=('1e-20','1e-25'),
        gates=dict(analytic_HP_relative=1e-25,cpp_value_relative=1e-12,double_FD_relative=1e-5),
        scope='32 frozen qualified BGN states at four previously tracked nodes; both carriers and both population partials. These are isolated constitutive partials, not a full assembled PhuMob Jacobian certification.',
        source='T-2022.03 equations 270-283 at 300 K, default arsenic; explicit chain rule and independent high precision centered differences.',
        acceptance_changes=False))
    a.write_csv(LOCAL/'input.csv',[{k:r[k] for k in ('id','nd','na','n','p')} for r in samples])
    d.matrix.freeze(OUT/'freeze.json',files+[OUT/'contract.json',LOCAL/'input.csv'])


def run():
    a.verify(OUT/'freeze.json')
    result=subprocess.run([str(p.LOCAL/'phumob_scalar_probe.exe'),str(LOCAL/'input.csv'),str(LOCAL/'cpp.csv'),'derivatives'],capture_output=True,text=True)
    a.write(OUT/'process.json',dict(exit_code=result.returncode,stdout=result.stdout,stderr=result.stderr));assert result.returncode==0
    cpp={x['id']:x for x in a.rows(LOCAL/'cpp.csv')};rows=[]
    with localcontext() as context:
        context.prec=100
        for sample in a.read(OUT/'contract.json')['samples']:
            state=[D.from_float(sample[k]) for k in ('nd','na','n','p')]
            for car in (0,1):
                mu,derivs,clamped=hp(state,car)
                value_error=float(abs(D(cpp[sample['id']]['mu_e' if car==0 else 'mu_h'])/mu-1))
                for pop in (0,1):
                    references=[]
                    for step in (D('1e-20'),D('1e-25')):
                        plus,minus=state.copy(),state.copy();plus[2+pop]*=step.exp();minus[2+pop]*=(-step).exp()
                        references.append((hp(plus,car)[0]-hp(minus,car)[0])/(2*step))
                    ref=references[-1];scale=abs(ref)
                    assert scale>0
                    hp_error=float(abs((derivs[pop]-ref)/ref))
                    step_error=float(abs((references[0]-ref)/ref))
                    for k,step in enumerate(STEPS):
                        computed=float(cpp[sample['id']][f'dlog_{car}_{pop}_{k}'])
                        error=float(abs((D.from_float(computed)-ref)/ref))
                        rows.append(dict(id=sample['id'],case=sample['case'],index=sample['index'],node=sample['node'],carrier=car,population=pop,
                            cross=(car!=pop),step_log_density=step,screening_clamped=clamped,cpp_value_relative=value_error,
                            analytic_hp_relative=hp_error,hp_step_relative=step_error,reference=float(ref),cpp_fd=computed,
                            fd_relative=error,relative_elasticity=float(abs(ref)/mu),fd_resolved=error<=1e-5))
    a.write_csv(OUT/'derivatives.csv',rows)
    summary=dict(samples=len(cpp),distinct_partials=len(rows)//len(STEPS),fd_checks=len(rows),
        cpp_value_max_relative=max(x['cpp_value_relative'] for x in rows),analytic_hp_max_relative=max(x['analytic_hp_relative'] for x in rows),hp_step_max_relative=max(x['hp_step_relative'] for x in rows),
        cpp_fd_failures=sum(not x['fd_resolved'] for x in rows),cpp_fd_zero_values=sum(x['cpp_fd']==0 for x in rows),
        unresolved_partials=sum(not any(r['fd_resolved'] for r in rows if (r['id'],r['carrier'],r['population'])==key) for key in {(r['id'],r['carrier'],r['population']) for r in rows}),
        native_values_not_yet_compared=True)
    a.write(OUT/'summary.json',summary)
    d.matrix.freeze(OUT/'evidence.json',[OUT/'freeze.json',LOCAL/'cpp.csv',OUT/'process.json',OUT/'derivatives.csv',OUT/'summary.json'])
    print(summary,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('prepare','run'))
    globals()[parser.parse_args().action]()
