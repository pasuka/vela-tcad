"""Independent fixed-state coefficient, port and changed-Jacobian preflight."""
import argparse
import math
from pathlib import Path
import numpy as np
import validate_simplemos_joint_geometry_20260907 as v
p=v.p;a=p.a;d=p.d;OUT=v.OUT;LOCAL=v.LOCAL


def jvp_metrics(path,key,arm,base=None):
    rows=a.rows(path);count=len(rows)//27;assert len(rows)==27*count
    grouped={}
    for r in rows:grouped.setdefault((int(r['direction']),int(r['block']),float(r['step_V'])),[]).append(r)
    results=[]
    for (mode,block,h),group in grouped.items():
        analytic=np.array([float(x['analytic']) for x in group]);fd=np.array([float(x['fd']) for x in group])
        error=float(np.linalg.norm(analytic-fd)/max(np.linalg.norm(analytic),1e-300));weak=(mode,block) in ((1,2),(2,1))
        identical='not_applicable'
        if weak and base is not None:
            ref=base[(mode,block,h)];identical=all(x['analytic']==y['analytic'] for x,y in zip(group,ref));assert identical
        passed=weak or h==1e-5 or error<=1e-4
        results.append(dict(key=key,arm=arm,input_block=mode,output_block=block,step_V=h,relative_error=error,
            unchanged_SRH_cross_block=weak,fixed_state_cross_block_bit_identical=identical,qualified=passed))
    return results,grouped


def main():
    a.verify(OUT/'freeze.json');contract=a.read(OUT/'contract.json');checks=[];jvp=[];ports=[];files=[Path(__file__).resolve(),OUT/'freeze.json']
    for c in contract['cases']:
        geo=d.matrix.spatial.m73.Geometry(c['device']);root=LOCAL/c['key']/'preflight';old=p.prior.LOCAL/'vela'/c['key']/'strict'
        baseR=d.ordered(root/'baseline/residual.csv',geo.count);oldR=d.ordered(old/'residual.csv',geo.count)
        baseE=a.rows(root/'baseline/edges.csv');oldE=a.rows(old/'edges.csv')
        assert all(x[k]==y[k] for x,y in zip(baseR,oldR) for k in ('psi_residual','phin_residual','phip_residual'))
        for x,y in zip(baseE,oldE):
            for k in ('couple_m','electron_flux','hole_flux','electron_mobility_m2_V_s','hole_mobility_m2_V_s'):assert x[k]==y[k]
        rawbase=jvp_metrics(root/'baseline/jvp.csv',c['key'],'baseline');jvp+=rawbase[0]
        eRat={int(x.split()[0]):float(x.split()[1]) for x in (Path(c['ratios'])/'edges.txt').read_text().splitlines()}
        baseline_abs=np.zeros((3,geo.count))
        for e in baseE:
            for block,car in [(1,'electron'),(2,'hole')]:
                i,j=int(e['node0']),int(e['node1']);f=abs(float(e[car+'_flux']));baseline_abs[block,i]+=f;baseline_abs[block,j]+=f
        for arm,(T,P) in v.AXES.items():
            r=d.ordered(root/arm/'residual.csv',geo.count);expected=d.ordered(root/arm/'expected.csv',geo.count)
            edges=a.rows(root/arm/'edges.csv');status=a.read(root/arm/'functional.status.json');assert status['exit_code']==0
            for i,(x,y) in enumerate(zip(edges,baseE)):
                ratio=(eRat[int(y['edge_id'])] if T else 1.)
                assert abs(float(x['couple_m'])-ratio*float(y['couple_m']))<=2e-14*max(abs(float(x['couple_m'])),abs(float(y['couple_m'])),1e-300)
                for car in ('electron','hole'):
                    assert x[car+'_mobility_m2_V_s']==y[car+'_mobility_m2_V_s']
                    f,oldf=float(x[car+'_flux']),float(y[car+'_flux']);assert abs(f-ratio*oldf)<=2e-13*max(abs(f),abs(oldf),1e-300)
            psi=np.array([float(x['psi_residual'])-float(y['psi_residual']) for x,y in zip(r,baseR)])
            expectedPsi=np.array([float(x['psi_source']) for x in expected])
            if P:perror=float(np.linalg.norm(psi-expectedPsi)/np.linalg.norm(expectedPsi));assert perror<=1e-8
            else:perror=0.;assert np.all(psi==0)
            maximum=0.
            for car,block,column in [('electron',1,'phin_residual'),('hole',2,'phip_residual')]:
                for i,(x,y,ex) in enumerate(zip(r,baseR,expected)):
                    actual=float(x[column])-float(y[column]);source=float(ex[car+'_source']);flow=float(ex[car+'_source_abs'])
                    bound=max(1e-8*flow+64*np.finfo(float).eps*baseline_abs[block,i],1e-25)
                    err=abs(actual-source);maximum=max(maximum,err/bound);assert err<=bound,(c['key'],arm,car,i,err,bound)
            current=status['current_A_per_um'];extractor=status['contact_current_extractor_A_per_um']
            porterror=abs(current/extractor-1);assert porterror<=1e-8
            prediction=next(x for x in a.rows(OUT/'predictions.csv') if x['key']==c['key'] and x['arm']==arm)
            oldcurrent=a.read(root/'baseline/functional.status.json')['current_A_per_um']
            assert abs(current-oldcurrent-float(prediction['unit_direct_A_per_um']))<=1e-8*max(abs(current),abs(oldcurrent))
            ports.append(dict(key=c['key'],arm=arm,functional_A_per_um=current,extractor_A_per_um=extractor,relative_error=porterror,qualified=True))
            checks.append(dict(key=c['key'],arm=arm,Poisson_source_relative_error=perror,carrier_source_error_over_bound=maximum,qualified=True))
            if arm!='baseline':jvp+=jvp_metrics(root/arm/'jvp.csv',c['key'],arm,rawbase[1])[0]
        files += [x for x in root.rglob('*') if x.is_file()]
        print('Fixed-state coefficient/port/JVP checked',c['key'],flush=True)
    a.write_csv(OUT/'preflight_checks.csv',checks);a.write_csv(OUT/'preflight_jvp.csv',jvp);a.write_csv(OUT/'preflight_ports.csv',ports)
    failed=[x for x in jvp if not x['qualified']]
    if failed:raise RuntimeError('Changed-block JVP preflight did not qualify: '+str(failed[:6]))
    files += [OUT/n for n in ('preflight_checks.csv','preflight_jvp.csv','preflight_ports.csv')]
    d.matrix.freeze(OUT/'preflight_evidence.json',files);print('All preflight gates qualified',flush=True)


if __name__=='__main__':main()
