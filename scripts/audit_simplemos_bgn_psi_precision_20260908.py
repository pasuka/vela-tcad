"""Independent Decimal-100 audit of actual BGN edges before a Jacobian repair.

The reference differentiates unfactorized density fluxes, not J or Jv output.
All numbers are per unit SG coefficient, so unit scaling cannot hide errors.
"""
import math
from decimal import Decimal as D, localcontext
from pathlib import Path
import simplemos_bgn_intrinsic_control_20260908 as prior

w=prior.w
a,d=w.a,w.d
LOCAL=w.REPO/'build-release/simplemos_bgn_psi_derivative_20260908'
OUT=w.REPO/'reference_tcad/simplemos_sentaurus2022/bgn_psi_derivative_20260908'
VT=1.380649e-23*300/1.602176634e-19


def bd(x):
    if abs(x)<1e-10:return 1.-x*.5+x*x/12
    if x>500:return x*math.exp(-x)
    if x < -500:return -x
    return x/math.expm1(x)


def dbd(x):
    if abs(x)<1e-8:return -.5+x/6-x*x*x/180
    if x>500:return (1-x)*math.exp(-x)
    if x < -500:return -1.
    e=math.expm1(x)
    return (e-x*(e+1))/(e*e)


def bp(x):
    return D(1) if x==0 else x/(x.exp()-1)


def dbp(x):
    if x==0:return D('-.5')
    e=x.exp()-1
    return (e-x*(e+1))/(e*e)


def density_flux(v,electron):
    n0,n1,p0,p1,q0,q1,vt=v
    eta=(p1-p0)/vt+(n1/n0 if electron else n0/n1).ln()
    left=n0*((p0-q0)/vt if electron else (q0-p0)/vt).exp()
    right=n1*((p1-q1)/vt if electron else (q1-p1)/vt).exp()
    return bp(-eta)*left-bp(eta)*right if electron else bp(eta)*left-bp(-eta)*right


def expanded(v,electron,high=False):
    ni0,ni1,psi0,psi1,q0,q1,vt=v
    ln=(lambda x:x.ln()) if high else math.log
    exp=(lambda x:x.exp()) if high else math.exp
    B,db=(bp,dbp) if high else (bd,dbd)
    eta=(psi1-psi0)/vt+ln(ni1/ni0 if electron else ni0/ni1)
    if electron:
        n0,n1=ni0*exp((psi0-q0)/vt),ni1*exp((psi1-q1)/vt)
        return [((db(-eta)+B(-eta))*n0+db(eta)*n1)/vt,
                (-db(-eta)*n0-(db(eta)+B(eta))*n1)/vt]
    p0,p1=ni0*exp((q0-psi0)/vt),ni1*exp((q1-psi1)/vt)
    return [(-(db(eta)+B(eta))*p0-db(-eta)*p1)/vt,
            (db(eta)*p0+(db(-eta)+B(-eta))*p1)/vt]


def factored(v,electron):
    ni0,ni1,psi0,psi1,q0,q1,vt=v
    eta=(psi1-psi0)/vt+math.log(ni1/ni0 if electron else ni0/ni1)
    if electron:
        f=bd(eta)*ni1*math.exp((psi1-q1)/vt)*math.expm1((q1-q0)/vt)
    else:
        f=bd(-eta)*ni1*math.exp((q1-psi1)/vt)*math.expm1((q0-q1)/vt)
    x=-eta if electron else eta
    log_derivative=dbd(x)/bd(x)
    return [f/vt*(1+log_derivative),-f/vt*log_derivative] if electron else [-f/vt*(1+log_derivative),f/vt*log_derivative]


def main():
    a.verify(prior.OUT/'comparison_evidence.json')
    jobs=[]
    for c in a.read(prior.OUT/'contract.json')['cases']:
        root=prior.LOCAL/'vela'/c['case']/'vela/vg_040'
        path=next(p/'acceptance_edges.csv' for p in sorted(root.glob('attempt_*')) if a.read(p/'result.json')['qualified'])
        jobs.append(dict(case=c['case'],path=str(path)))
    a.write(OUT/'precision_contract.json',dict(jobs=jobs,precision=100,finite_difference_step_V='1e-25',
        selection='Per case and carrier choose 12 silicon edges with largest absolute expanded-versus-factored double derivative difference. Selection is deterministic; no result threshold selects edges.',
        reference='Central differences of independently reconstructed unfactorized SG density flux at Decimal-100; also evaluate original expanded derivative formula at Decimal-100.',
        parameters='Exact binary64 values of exported ni/psi/quasi-Fermi potentials, promoted to Decimal; coefficient=1. Physical-potential CSV audit is not a replay of internally split quasi-Fermi references.',
        gates=dict(expanded_high_precision_relative=1e-30,factored_double_relative=1e-7),
        purpose='Distinguish algebraic omissions from subtraction loss. Full production Jv is independently required.'))
    d.matrix.freeze(OUT/'precision_freeze.json',[Path(__file__).resolve(),OUT/'precision_contract.json',prior.OUT/'comparison_evidence.json']+[Path(j['path']) for j in jobs])
    rows=[]
    for job in jobs:
        edges=a.rows(Path(job['path']))
        for electron in (True,False):
            candidates=[]
            for r in edges:
                if float(r['ni0_m3'])<=0 or float(r['ni1_m3'])<=0 or float(r['couple_m'])==0:continue
                q='phin' if electron else 'phip'
                values=[float(r[k]) for k in ('ni0_m3','ni1_m3','psi0_V','psi1_V',q+'0_V',q+'1_V')]+[VT]
                if any(abs((values[i+2]-values[i+4])/VT)>=500 for i in (0,1)):continue
                old,new=expanded(values,electron),factored(values,electron)
                candidates.append((max(abs(x-y) for x,y in zip(old,new)),int(r['edge_id']),r,values,old,new))
            for _,edge,r,values,old,new in sorted(candidates,key=lambda c:(-c[0],c[1]))[:12]:
                with localcontext() as context:
                    context.prec=100
                    hp=[D.from_float(x) for x in values]
                    h=D('1e-25')
                    hp_expanded=expanded(hp,electron,True)
                    for end in (0,1):
                        plus,minus=hp.copy(),hp.copy()
                        plus[2+end]+=h;minus[2+end]-=h
                        reference=(density_flux(plus,electron)-density_flux(minus,electron))/(2*h)
                        # Exact flat-qF equilibrium derivatives are zero; residual
                        # Decimal cancellation is recorded as an absolute error.
                        flat=values[4]==values[5]
                        scale=abs(reference) if not flat else D(1)
                        old_error=float(abs(D.from_float(old[end])-reference)/scale)
                        new_error=float(abs(D.from_float(new[end])-reference)/scale)
                        high_error=float(abs(hp_expanded[end]-reference)/scale)
                        rows.append(dict(case=job['case'],carrier='electron' if electron else 'hole',edge_id=edge,
                            node0=int(r['node0']),node1=int(r['node1']),endpoint=end,flat_qf=flat,
                            reference=float(reference),expanded_double=old[end],factored_double=new[end],
                            expanded_double_error=old_error,factored_double_error=new_error,expanded_hp_error=high_error,
                            error_metric='absolute at flat qF' if flat else 'relative to independent HP finite difference'))
    w.v.csv_union(OUT/'precision_edges.csv',rows)
    summary=dict(checks=len(rows),expanded_hp_max_error=max(r['expanded_hp_error'] for r in rows),
        expanded_double_max_error=max(r['expanded_double_error'] for r in rows),factored_double_max_error=max(r['factored_double_error'] for r in rows),
        high_precision_formula_pass=all(r['expanded_hp_error']<=1e-30 for r in rows),
        factored_pass=all(r['factored_double_error']<=1e-7 for r in rows))
    a.write(OUT/'precision_summary.json',summary)
    d.matrix.freeze(OUT/'precision_evidence.json',[OUT/'precision_freeze.json',OUT/'precision_edges.csv',OUT/'precision_summary.json'])
    print(summary,flush=True)


if __name__=='__main__':main()
