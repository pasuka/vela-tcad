"""Independent high-precision SG/SRH accounting on frozen binary inputs."""
from decimal import Decimal,localcontext,getcontext
from pathlib import Path
import math
import numpy as np
import audit_simplemos_minority_residual_20260906 as v

a=v.a;d=v.d;D=Decimal


def exact(x):return D.from_float(float(x))


def em1(x):
    if abs(x)>=D('.01'):return x.exp()-1
    term=x;total=x;k=2
    while True:
        term=term*x/k;new=total+term
        if new==total:return total
        total=new;k+=1


def bern(x):return D(1) if x==0 else x/em1(x)
def clamp(x):return max(D(-500),min(D(500),x))


def edge_flux(e,nodes,carrier,mode):
    c=e['ncoef' if carrier=='electron' else 'pcoef']
    if c==0:return D(0)
    i=int(e['node0']);j=int(e['node1']);vt=e['Vt'];ni0=e['ni0'];ni1=e['ni1']
    prefix='e' if carrier=='electron' else 'h'
    if mode=='kernel':psi0=e[prefix+'psi0'];psi1=e[prefix+'psi1'];q0=e[prefix+'phi0'];q1=e[prefix+'phi1']
    else:
        p=nodes[i];q=nodes[j];psi0=p['psi'];psi1=q['psi'];q0=p[prefix+'ref']+p[prefix+'inc'];q1=q[prefix+'ref']+q[prefix+'inc']
    electron=carrier=='electron';eta=(psi1-psi0)/vt
    z0=(psi0-q0)/vt if electron else (q0-psi0)/vt
    z1=(psi1-q1)/vt if electron else (q1-psi1)/vt
    b0=clamp(z0);b1=clamp(z1)
    if ni0!=ni1:
        # Current no-BGN production policy uses density SG at a material interface.
        left=bern(-eta if electron else eta)*ni0*b0.exp()
        right=bern(eta if electron else -eta)*ni1*b1.exp()
        return c*(left-right)/e['scale']
    if ni0<=0:return D(0)
    drive=((q1-q0) if electron else (q0-q1))/vt+(b0-z0)-(b1-z1)
    return c*bern(eta if electron else -eta)*ni1*b1.exp()*em1(drive)/e['scale']


def srh(node,mode):
    ni=node['ni'];vt=node['Vt']
    if ni<=0:return D(0)
    if mode=='kernel':n=node['nsrh'];p=node['p'];dp=node['dphi']
    else:
        phin=node['eref']+node['einc'];phip=node['href']+node['hinc'];dp=phip-phin
        n=ni*clamp((node['psi']-phin)/vt).exp();p=ni*clamp((phip-node['psi'])/vt).exp()
    excess=ni*ni*em1(clamp(dp/vt));den=node['taup']*(n+ni)+node['taun']*(p+ni)
    assert abs(excess)<D('1e60') and den>D('1e-100'),'Outside the frozen unclipped SRH reference branch'
    rate=excess/den
    return rate*node['volume']*node['source_factor']/node['scale']


def evaluate(edges,nodes,precision):
    with localcontext() as ctx:
        ctx.prec=precision;result={}
        for mode in ('kernel','partition'):
            source=[srh(n,mode) for n in nodes]
            for carrier,field in (('electron','nflux'),('hole','pflux')):
                flux=[D(0) for n in nodes];absolute=[D(0) for n in nodes];ef=[]
                for e in edges:
                    f=edge_flux(e,nodes,carrier,mode);ef.append(f);i=int(e['node0']);j=int(e['node1'])
                    flux[i]+=f;flux[j]-=f;absolute[i]+=abs(f);absolute[j]+=abs(f)
                result[mode,carrier]=dict(flux=flux,source=source,absolute=absolute,edges=ef)
        return result


def run():
    getcontext().prec=120
    a.verify(v.OUT/'freeze.json');a.verify(v.m.OUT/'validation_evidence.json')
    rows=[];summary=[];edge_rows=[];files=[Path(__file__).resolve(),v.OUT/'freeze.json']
    for job in a.read(v.OUT/'contract.json')['probes']:
        path=Path(job['config']).parent;src=Path(job['source']);cfg=a.read(path/'config.json')
        assert cfg['solver']['bandgap_narrowing']['model']=='none' and cfg['solver']['recombination']==['srh']
        nodes=[{k:exact(x) for k,x in r.items()} for r in a.rows(path/'kernel_nodes.csv')]
        edges=[{k:exact(x) for k,x in r.items()} for r in a.rows(path/'kernel_edges.csv')]
        terms=d.ordered(src/'all_row.csv',len(nodes));geo=d.matrix.spatial.m73.Geometry(job['device'])
        masks=d.matrix.spatial.old.m78.supports(job['device'],geo,.05)[0];mask=masks['all_si'].copy();mask[geo.contact_nodes]=False
        results={p:evaluate(edges,nodes,p) for p in (60,100)};group=[]
        for carrier,field in (('electron','nflux'),('hole','pflux')):
            prod=[D(0) for n in nodes]
            for e in edges:prod[int(e['node0'])]+=e[field];prod[int(e['node1'])]-=e[field]
            for i in np.flatnonzero(mask):
                t=terms[i];f=float(t[carrier+'_flux']);source=float(t[carrier+'_recombination']);scale=max(float(t[carrier+'_flux_abs_sum']),abs(source))
                assert scale>0
                recon=abs(float(prod[i])+source-float(t[carrier+'_residual']))/scale
                assert recon<=1e-10,(job['tag'],i,carrier,recon)
                row=dict(tag=job['tag'],case=job['case'],vg=job['vg'],carrier=carrier,node=int(i),channel=bool(masks['channel'][i]),
                    density_cm3=float(nodes[i]['n' if carrier=='electron' else 'p']),baseline_ratio=abs(float(t[carrier+'_residual']))/scale,
                    baseline_scale=scale,assembly_reconstruction_scaled_error=recon)
                for mode in ('kernel','partition'):
                    q=results[100][mode,carrier];fhi=q['flux'][i];shi=q['source'][i];rhi=fhi+shi;s=max(q['absolute'][i],abs(shi))
                    low=results[60][mode,carrier];agreement=float(abs((low['flux'][i]+low['source'][i])-rhi))/scale
                    assert agreement<=1e-20,(job['tag'],i,carrier,mode,agreement)
                    row.update({mode+'_flux':float(fhi),mode+'_source':float(shi),mode+'_residual':float(rhi),mode+'_ratio':float(abs(rhi)/s) if s else None,
                        mode+'_zero_scale':s==0,mode+'_residual_change_over_baseline_scale':float(rhi-exact(t[carrier+'_residual']))/scale,
                        mode+'_precision_agreement':agreement})
                group.append(row);rows.append(row)
        worst=max((r for r in group if r['carrier']=='hole'),key=lambda r:r['baseline_ratio']);wi=worst['node']
        for k,e in enumerate(edges):
            if wi not in (int(e['node0']),int(e['node1'])):continue
            edge_rows.append(dict(tag=job['tag'],worst_hole_node=wi,edge=int(e['edge']),node0=int(e['node0']),node1=int(e['node1']),
                production_flux=float(e['pflux']),kernel_flux=float(results[100]['kernel','hole']['edges'][k]),partition_flux=float(results[100]['partition','hole']['edges'][k]),
                kernel_qf_difference=float(e['hphi0']-e['hphi1']),partition_qf_difference=float(nodes[int(e['node0'])]['href']+nodes[int(e['node0'])]['hinc']-nodes[int(e['node1'])]['href']-nodes[int(e['node1'])]['hinc'])))
        s=dict(tag=job['tag'],case=job['case'],vg=job['vg'],initialization=job['initialization'],rows=len(group),worst_hole_node=wi,
            baseline_violations=sum(r['baseline_ratio']>1e-6 for r in group),baseline_max=max(r['baseline_ratio'] for r in group))
        for mode in ('kernel','partition'):
            s[mode+'_violations']=sum(r[mode+'_ratio'] is None or r[mode+'_ratio']>1e-6 for r in group)
            s[mode+'_max']=max(r[mode+'_ratio'] or 0 for r in group)
            s[mode+'_max_residual_change_scaled']=max(abs(r[mode+'_residual_change_over_baseline_scale']) for r in group)
        summary.append(s);files += [path/'kernel_nodes.csv',path/'kernel_edges.csv',path/'config.status.json',src/'all_row.csv']
        print(job['tag'],s['baseline_violations'],s['kernel_violations'],s['partition_violations'],flush=True)
    a.write_csv(v.OUT/'precision_rows.csv',rows);a.write_csv(v.OUT/'precision_summary.csv',summary);a.write_csv(v.OUT/'worst_hole_edges.csv',edge_rows)
    a.write(v.OUT/'precision_result.json',dict(input_hashes={a.rel(p):a.sha(p) for p in sorted(set(files))},states=len(summary),rows=len(rows),
        max_reconstruction=max(r['assembly_reconstruction_scaled_error'] for r in rows),max_precision_agreement=max(r[x+'_precision_agreement'] for r in rows for x in ('kernel','partition')),
        boundary='Fixed-state arithmetic diagnosis; high-precision residual evaluation is not a high-precision self-consistent solve. Mobility/lifetimes remain at production double values.',production_changes=False))


if __name__=='__main__':run()
