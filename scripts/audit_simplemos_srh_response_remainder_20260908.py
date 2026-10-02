"""Read-only attribution of the qualified n23/low-Vd response discrepancy."""
import argparse
from pathlib import Path
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import splu
import calibrate_simplemos_srh_hotspots_20260908 as s

a=s.a;d=s.d;t=s.t;v=s.v;p=s.p
OUT=s.OUT/'qualified_remainder';LOCAL=s.LOCAL/'qualified_remainder'


def prepare():
    a.verify(s.OUT/'evidence.json')
    case=next(c for c in a.read(s.OUT/'contract.json')['cases'] if c['device']=='n23' and c['vd']==.05)
    dc=[r for r in a.rows(s.OUT/'dc.csv') if r['key']==case['key']];assert len(dc)==5 and all(r['qualified']=='True' for r in dc)
    files=[Path(__file__).resolve(),s.OUT/'evidence.json'];jobs=[]
    for job in case['jobs']:
        if job['label']=='zero':continue
        original=Path(job['dest']);dest=LOCAL/job['label'];cfg=a.read(original/'config.json')
        files+=v.probes(cfg,dest,original/'state.csv');files.append(original/'state.csv')
        jobs.append(dict(label=job['label'],alpha=job['alpha'],dest=str(dest),original=str(original)))
    a.write(OUT/'contract.json',dict(case=case,jobs=jobs,scope='Read-only residuals at four already-qualified frozen perturbed states; no new DC or acceptance changes.',
        identity='J*(odd state - alpha*tangent) versus residual odd part; split Poisson/electron/hole residual projections. This is local linear attribution, not a new calibration pass.'))
    d.matrix.freeze(OUT/'freeze.json',files+[OUT/'contract.json'])


def run():
    a.verify(OUT/'freeze.json');contract=a.read(OUT/'contract.json');case=contract['case']
    for job in contract['jobs']:
        status=s.execute(Path(job['dest'])/'functional.json',case,job['alpha']);assert status['exit_code']==0,status
    geo,mask=v.m.previous.prior.support(case);N=geo.count;ids=np.flatnonzero(mask);fixed=s.LOCAL/case['key']/'fixed/zero'
    J=t.matrix(fixed/'jacobian.csv',N);wide=J.tocsr().astype(np.longdouble)
    scale=1/np.maximum(np.asarray(abs(J).max(axis=1).toarray()).ravel(),1e-300);lu=splu((diags(scale)@J).tocsc())
    V=a.read(fixed/'jacobian.status.json')['potential_scale_V']
    def solve(rhs):
        dx=lu.solve(scale*rhs)
        for _ in range(4):dx-=lu.solve(scale*np.asarray(wide@dx.astype(np.longdouble)-rhs.astype(np.longdouble),dtype=float))
        linear=np.linalg.norm(np.asarray(wide@dx.astype(np.longdouble)-rhs.astype(np.longdouble),dtype=float))/max(np.linalg.norm(rhs),1e-300)
        return dx,linear
    def field(dx):return (dx.reshape(3,N)[:,ids]*V)
    base=d.ordered(Path(case['baseline'])/'state.csv',N);terms=d.ordered(fixed/'terms.csv',N);source=np.zeros(3*N)
    for site in case['hotspots']:
        for block,carrier in ((1,'electron'),(2,'hole')):source[block*N+site['node']]=float(terms[site['node']][carrier+'_recombination'])*(site['ratio']-1.)
    tangent,lin=solve(-source);rows=[];parts=[];fields=[]
    for amp,alpha in [('full',.001),('half',.0005)]:
        plus=next(j for j in contract['jobs'] if j['label']=='plus_'+amp);minus=next(j for j in contract['jobs'] if j['label']=='minus_'+amp)
        delta=(t.physical_delta(d.ordered(Path(plus['original'])/'state.csv',N),base,np.arange(N))-t.physical_delta(d.ordered(Path(minus['original'])/'state.csv',N),base,np.arange(N)))/2
        odd=(s.c.old.R(Path(plus['dest']),N)-s.c.old.R(Path(minus['dest']),N))/2
        error=delta[:,ids]-alpha*field(tangent);pred=alpha*field(tangent);rproject,projlin=solve(odd);projection=field(rproject)
        unaccounted=np.linalg.norm(error-projection)/np.linalg.norm(pred)
        rows.append(dict(amplitude=amp,original_prediction_relative=float(np.linalg.norm(error)/np.linalg.norm(pred)),residual_projected_relative=float(np.linalg.norm(projection)/np.linalg.norm(pred)),
            remaining_relative=float(unaccounted),projection_fraction=float(np.sum(error*projection)/np.sum(error*error)),
            source_linear_relative=float(lin),projection_linear_relative=float(projlin),calibration_requalified=False))
        for block,name in enumerate(('poisson','electron','hole')):
            rhs=np.zeros(3*N);rhs[block*N:(block+1)*N]=odd[block*N:(block+1)*N];part,partlin=solve(rhs);pf=field(part)
            parts.append(dict(amplitude=amp,block=name,projection_relative=float(np.linalg.norm(pf)/np.linalg.norm(pred)),error_direction_fraction=float(np.sum(error*pf)/np.sum(error*error)),linear_relative=float(partlin)))
        for j,node in enumerate(ids):
            row=dict(amplitude=amp,node=int(node))
            for block,name in enumerate(('psi','phin','phip')):row.update({name+'_error_V':error[block,j],name+'_projection_V':projection[block,j],name+'_remaining_V':error[block,j]-projection[block,j]})
            fields.append(row)
    for name,data in [('metrics',rows),('blocks',parts),('fields',fields)]:a.write_csv(OUT/(name+'.csv'),data)
    d.matrix.freeze(OUT/'evidence.json',[OUT/'freeze.json',OUT/'metrics.csv',OUT/'blocks.csv',OUT/'fields.csv']+[f for f in LOCAL.rglob('*') if f.is_file()])
    print(rows,flush=True);print(parts,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('prepare','run'));globals()[parser.parse_args().action]()
