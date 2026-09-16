"""Decimal audit of recorded hole-row operands; no nonlinear updates or gate changes.

Precision applies to arithmetic on the EXACT captured binary64 operands. Material
statistics and mobility remain frozen, so this is not a full arbitrary-precision
Fermi/IALMob residual. Decimal 70/100 digit agreement checks truncation error.
"""
import argparse
import json
import math
from decimal import Decimal, localcontext
from pathlib import Path


def D(x):
    return Decimal.from_float(x) if isinstance(x,float) else Decimal(x)


def sg(edge, row, rounded_qf=False, frozen_g=False):
    a,b=edge['sa'],edge['sb'];pa,pb=edge['pa'],edge['pb']
    vt=D(row['kb'])/D(row['q'])*(D(a['T'])+D(b['T']))/2
    lf=(D(pb['p'])/D(pb['Nv'])).ln()-(D(pa['p'])/D(pa['Nv'])).ln()
    # Keep the production near-equal-eta fallback as a frozen operand.
    g=D(edge['g']) if frozen_g or abs(edge['logF'])<=1e-8 else (D(pb['ep'])-D(pa['ep']))/lf
    if g<=0:g=D(edge['g'])
    dq=D(b['rp'])-D(a['rp'])+D(b['fp'])-D(a['fp'])
    if rounded_qf:dq=D((b['rp']-a['rp'])+b['fp']-a['fp'])
    y=-dq/(vt*g);arg=(D(pb['p'])/D(pa['p'])).ln()+y
    bernoulli=arg/(arg.exp()-1) if arg else D(1)
    return D(row['q'])*D(edge['mu'])*D(edge['weight'])*vt*g*D(pb['p'])*bernoulli*(y.exp()-1)


def source(row, with_generation):
    s=row['state'];p=row['properties'];n=D(p['n']);h=D(p['p']);ni=D(p['ni'])
    vt=D(row['kb'])/D(row['q'])*D(s['T'])
    split=(D(s['rp'])-D(s['rn'])+D(s['fp'])-D(s['fn']))/vt
    gn=((n/D(p['Nc'])).ln()-D(p['en'])).exp()
    gp=((h/D(p['Nv'])).ln()-D(p['ep'])).exp()
    excess=n*h*(1-(-split).exp()) if split>=0 else gn*gp*ni*ni*(split.exp()-1)
    srh=excess/(D(p['tp'])*(n+gn*ni)+D(p['tn'])*(h+gp*ni))
    auger=(D(p['Cn'])*n+D(p['Cp'])*h)*excess if with_generation or split>0 else D(0)
    return D(row['q_area'])*(srh+auger)


def contact(row, potential_shift=Decimal(0)):
    if 'contact' not in row:return D(0)
    c=row['contact'];s=row['state']
    delta=(D(s['rp'])-D(c['bias'])+D(s['fp'])-(D(s['psi'])-D(c['neutral_potential'])-potential_shift))/D(c['vt'])
    if abs(delta)>=D(1e-6):raise ValueError('Contact audit requires the recorded small-delta branch')
    return D(c['coefficient'])*D(c['Nv'])*delta*(D(c['df'])+delta*D(c['ddf'])/2)


def audit(row, config, with_generation, precision):
    with localcontext() as ctx:
        ctx.prec=precision
        edges=row['edges'];signs=[1 if e['a']==row['node'] else -1 for e in edges]
        currents=[D(e['current']) for e in edges]
        captured_source=D(row['source']);scale=D(row['current_scale'])
        floor=D(config.get('scale_floor',0))*scale
        captured_contact=D(row.get('contact',{}).get('outward',0))
        def metric(values, src, boundary=captured_contact):
            residual=sum((s*v for s,v in zip(signs,values)),src)+boundary
            absolute=sum((abs(v) for v in values),abs(boundary))
            denominator=max(absolute,abs(src),floor)
            return dict(residual=float(residual),flux_abs=float(absolute),denominator=float(denominator),
                ratio=float(abs(residual)/denominator) if denominator else 0.)
        summation=metric(currents,captured_source)
        sr=D(row['q_area'])*(D(row['properties']['srh'])+D(row['properties']['auger']))
        precise=[sg(e,row) for e in edges]
        rounded=[sg(e,row,rounded_qf=True) for e in edges]
        frozen=[sg(e,row,frozen_g=True) for e in edges]
        source_precise=source(row,with_generation)
        contact_precise=contact(row)
        neutral=row.get('contact',{}).get('neutral_potential',0.)
        neutral_ulp=D(2)**(math.frexp(abs(neutral))[1]-53 if neutral else -1074)
        production_denominator=max(abs(row['flux_abs']/row['current_scale']),abs(row['source']/row['current_scale']),config.get('scale_floor',0))
        production_ratio=abs(row['residual']/row['current_scale'])/production_denominator
        # Accumulation round-off must not be confused with edge evaluation error.
        details=[]
        for e,v,rounded_v in zip(edges,precise,rounded):
            dq=D(e['sb']['rp'])-D(e['sa']['rp'])+D(e['sb']['fp'])-D(e['sa']['fp'])
            details.append(dict(id=e['id'],a=e['a'],b=e['b'],current=e['current'],decimal_current=float(v),
                rounded_qf_current=float(rounded_v),exact_qf_difference=str(dq),
                production_qf_difference=(e['sb']['rp']-e['sa']['rp'])+e['sb']['fp']-e['sa']['fp'],
                absolute_error=float(abs(D(e['current'])-v))))
        return dict(node=row['node'],precision=precision,production_residual=row['residual'],production_ratio=production_ratio,
            captured_sum=summation,source_rate_sum=metric(currents,sr),sg_rounded_qf=metric(rounded,captured_source),
            sg_frozen_g=metric(frozen,captured_source),sg_recomputed=metric(precise,captured_source),
            recombination_recomputed=metric(currents,source_precise),all_recomputed=metric(precise,source_precise),
            contact_recomputed=metric(currents,captured_source,contact_precise),
            all_with_contact=metric(precise,source_precise,contact_precise),
            neutral_plus_ulp=metric(currents,captured_source,contact(row,neutral_ulp)),
            neutral_minus_ulp=metric(currents,captured_source,contact(row,-neutral_ulp)),neutral_ulp=float(neutral_ulp),
            source=row['source'],decimal_source=float(source_precise),contact=row.get('contact'),state=row['state'],edges=details)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    args=p.parse_args();results=[]
    for path in sorted(args.root.glob('audit_*/output.json')):
        data=json.loads(path.read_text());cfg=json.loads(path.with_name('input.json').read_text())
        for row in data['hole_row_audit']:
            lo=audit(row,cfg['electrical_gate_solver']['carrier_row_convergence'],data['auger_with_generation'],70)
            hi=audit(row,cfg['electrical_gate_solver']['carrier_row_convergence'],data['auger_with_generation'],100)
            for key in ('captured_sum','source_rate_sum','sg_rounded_qf','sg_frozen_g','sg_recomputed','recombination_recomputed','all_recomputed','contact_recomputed','all_with_contact','neutral_plus_ulp','neutral_minus_ulp'):
                if lo[key]!=hi[key]:raise ValueError('Decimal precision did not stabilize: '+key)
            # Directly reconstruct the exact production accumulation order.
            ordered=row['source'];absolute=0.
            for e in row['edges']:
                ordered+=(1 if e['a']==row['node'] else -1)*e['current'];absolute+=abs(e['current'])
            if 'contact' in row:
                ordered+=row['contact']['outward'];absolute+=abs(row['contact']['outward'])
            if ordered!=row['residual'] or absolute!=row['flux_abs']:raise ValueError('Incomplete row terms')
            hi.update(case=path.parent.name,precision_stable=True,production_order_exact=True)
            results.append(hi)
    if not results:raise ValueError('No audited states')
    (args.root/'roundoff_findings.json').write_text(json.dumps(results,indent=2,allow_nan=False),encoding='utf-8')
    for r in results:print(r['case'],r['node'],r['production_ratio'],r['captured_sum']['ratio'],r['sg_recomputed']['ratio'],r['all_with_contact']['ratio'],r['neutral_plus_ulp']['ratio'])


if __name__=='__main__':main()
