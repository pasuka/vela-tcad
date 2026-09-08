"""Separate four-site SRH direction: legacy solver, frozen joint baselines."""
import argparse
import copy
import subprocess
import math
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import splu
import build_simplemos_stable_merit_20260908 as num

t=num.t;c=t.c;a=t.a;d=t.d;v=t.v;p=t.p
LOCAL=p.REPO/'build-release/simplemos_srh_hotspots_20260908'
OUT=p.REPO/'reference_tcad/simplemos_sentaurus2022/srh_hotspots_20260908'
RUNNER=LOCAL/'runner.exe'


def build():
    a.verify(num.audit.OUT/'final_evidence.json');a.verify(t.OUT/'build_evidence.json')
    LOCAL.mkdir(parents=True,exist_ok=False);OUT.mkdir(parents=True,exist_ok=False)
    source=(t.LOCAL/'CoupledDDAssembler.cpp').read_text(encoding='utf-8');assert source.count(t.b.HOOK)==1
    extra=t.b.HOOK.replace('VELA_CANDIDATE_SRH_','VELA_NEIGHBORHOOD_SRH_').replace('count==0||count>2','count!=4')
    target=LOCAL/'CoupledDDAssembler.cpp';target.write_text(source.replace(t.b.HOOK,t.b.HOOK+'\n'+extra),encoding='utf-8')
    cmd=list(a.read(t.LOCAL/'compile_command.json'));old=cmd[cmd.index('-o')+1];cmd[cmd.index('-c')+1]=str(target);cmd[cmd.index('-o')+1]=str(LOCAL/'CoupledDDAssembler.o')
    a.write(LOCAL/'compile_command.json',cmd);r=subprocess.run(cmd,cwd=p.REPO/'build-release',env=v.environment(),capture_output=True,encoding='utf-8',errors='replace');(LOCAL/'build.log').write_text(r.stdout+r.stderr,encoding='utf-8');assert r.returncode==0,r.stderr[-3000:]
    link=[str(LOCAL/'CoupledDDAssembler.o') if x==old else x for x in a.read(t.LOCAL/'link_command.json')];link[link.index('-o')+1]=str(RUNNER)
    a.write(LOCAL/'link_command.json',link);r=subprocess.run(link,env=v.environment(),capture_output=True,encoding='utf-8',errors='replace');(LOCAL/'link.log').write_text(r.stdout+r.stderr,encoding='utf-8');assert r.returncode==0,r.stderr[-3000:]
    d.matrix.freeze(OUT/'build_evidence.json',[Path(__file__).resolve(),num.audit.OUT/'final_evidence.json',t.OUT/'build_evidence.json']+[x for x in LOCAL.rglob('*') if x.is_file()]+[Path(x) for x in link if x.endswith('.o')]);print('Built four-site SRH-only overlay with original Newton comparison',flush=True)


def env(case,alpha):
    e=t.env(case,'joint',1.);e.update(VELA_NEIGHBORHOOD_SRH_VOLUMES=str(LOCAL/(case['device']+'_sites.txt')),VELA_NEIGHBORHOOD_SRH_ALPHA=format(alpha,'.17g'))
    return e


def execute(path,case,alpha):return v.execute(path,RUNNER,env(case,alpha))


def prepare():
    a.verify(OUT/'build_evidence.json');cases=a.read(t.OUT/'contract.json')['cases'];geo23=c.s.prior.old.l.g.geometry('n23')[0];scope=[];files=[Path(__file__).resolve(),OUT/'build_evidence.json',RUNNER]
    maps={}
    for dev in ('n19','n23'):
        geo,xy,cells,vol,*_=c.s.prior.old.l.g.geometry(dev);sites=[]
        for tag in (794,795,1056,1087):
            target=np.array(geo23.coords[tag]);dist=[np.linalg.norm(np.array(geo.coords[i])-target) for i in range(geo.count)];node=int(np.argmin(dist));assert dist[node]<=1e-8 and node in geo.free
            ratio=float(vol['Si'][node]/geo.volumes['all_cell'][node]);assert ratio>0
            sites.append(dict(node=node,ratio=ratio));scope.append(dict(device=dev,tag=tag,node=node,distance_um=dist[node],x_um=geo.coords[node][0],y_um=geo.coords[node][1],ratio=ratio))
        maps[dev]=sites;path=LOCAL/(dev+'_sites.txt');path.write_text(f'{geo.count} 4\n'+''.join(f"{x['node']} {x['ratio']:.17g}\n" for x in sites));files.append(path)
    for case in cases:
        case['hotspots']=maps[case['device']];assert not ({x['node'] for x in case['hotspots']}&{x['node'] for x in case['mapped_nodes'].values()})
        row=next(x for x in a.rows(num.audit.prior.OUT/'selected_states.csv') if x['key']==case['key'] and x['axis']=='joint' and x['label']=='finite');assert row['qualified']=='True'
        base=p.REPO/row['path'];case['baseline']=str(base);files.append(base/'state.csv')
        cfg=a.read(base/'config.json');cfg['solver'].pop('local_update_diagnostics',None)
        case['jobs']=[]
        for label,alpha in [('zero',0.),('plus_full',.001),('minus_full',-.001),('plus_half',.0005),('minus_half',-.0005)]:
            dest=LOCAL/case['key']/label;deck=copy.deepcopy(cfg);deck.update(state_file=str(base/'state.csv'),output_state_file=str(dest/'state.csv'))
            a.write(dest/'config.json',deck);v.post_config(deck,dest);files += [dest/(n+'.json') for n in ('config','all_row','acceptance_edges')]
            case['jobs'].append(dict(dest=str(dest),label=label,alpha=alpha))
        for label in ('zero','finite_fixed'):
            dest=LOCAL/case['key']/'fixed'/label;files+=v.probes(cfg,dest,base/'state.csv')
            for name,typ in [('jacobian','parameter_jacobian'),('jvp','joint_geometry_jvp')]:
                deck=a.read(dest/'functional.json');deck.update(simulation_type=typ,output_csv=str(dest/(name+'.csv')));a.write(dest/(name+'.json'),deck);files.append(dest/(name+'.json'))
    a.write_csv(OUT/'scope.csv',scope);a.write(OUT/'contract.json',dict(cases=cases,DC=20,source='Four spatially matched previously audited hotspot sites; interpolate each current SRH volume toward independent signed Si volume using one common alpha. Existing two joint-volume anchors stay fixed.',
        solver='Original production residual/merit/state algorithms; stable-merit candidate NOT enabled. Both NWell and Vd; no new self-consistent finite replacement.',
        response='Full physical psi/phin/phip vector on all 907 free Si nodes; independent complete J tangent versus plus/minus .001 and .0005.',
        gates=dict(row=1e-6,kcl=1e-8,source_relative=1e-5,linear_relative=1e-8,weak_J_relative=1e-12,prediction=1e-3,two_amplitude=1e-3,even_over_odd=.01,signal_over_drift=100),native_source_response_calibrated=False))
    d.matrix.freeze(OUT/'freeze.json',files+[OUT/'scope.csv',OUT/'contract.json']);print('Frozen 4-site direction / 20 DC',flush=True)


def preflight():
    a.verify(OUT/'freeze.json');checks=[];tangents=[]
    for case in a.read(OUT/'contract.json')['cases']:
        geo,mask=v.m.previous.prior.support(case);N=geo.count;root=LOCAL/case['key']/'fixed'
        for label,alpha in [('zero',0.),('finite_fixed',1.)]:
            for name in ('functional','edges','terms','jacobian','jvp'):
                status=execute(root/label/(name+'.json'),case,alpha);assert status['exit_code']==0,status
        zero=root/'zero';changed=root/'finite_fixed';R=c.old.R(zero,N);original=c.old.R(Path(case['baseline'])/'post',N);terms=d.ordered(zero/'terms.csv',N)
        source=np.zeros(3*N)
        for site in case['hotspots']:
            for block,car in ((1,'electron'),(2,'hole')):source[block*N+site['node']]=float(terms[site['node']][car+'_recombination'])*(site['ratio']-1.)
        source_error=float(np.linalg.norm(c.old.R(changed,N)-R-source)/np.linalg.norm(source))
        J=t.matrix(zero/'jacobian.csv',N);J1=t.matrix(changed/'jacobian.csv',N);weak=0.
        for ob,ib in ((1,2),(2,1)):
            old=J[ob*N:(ob+1)*N,ib*N:(ib+1)*N].toarray();new=J1[ob*N:(ob+1)*N,ib*N:(ib+1)*N].toarray();expected=old.copy()
            for site in case['hotspots']:expected[site['node'],:]*=site['ratio']
            nz=expected!=0;assert np.all(new[~nz]==0);weak=max(weak,float(np.max(abs(new[nz]/expected[nz]-1))))
        scale=1/np.maximum(np.asarray(abs(J).max(axis=1).toarray()).ravel(),1e-300);lu=splu((diags(scale)@J).tocsc());wide=J.tocsr().astype(np.longdouble);dx=lu.solve(-scale*source)
        for _ in range(4):dx-=lu.solve(scale*np.asarray(wide@dx.astype(np.longdouble)+source.astype(np.longdouble),dtype=float))
        linear=float(np.linalg.norm(np.asarray(wide@dx.astype(np.longdouble)+source.astype(np.longdouble),dtype=float))/np.linalg.norm(source));V=a.read(zero/'jacobian.status.json')['potential_scale_V']
        for i in np.flatnonzero(mask):tangents.append(dict(key=case['key'],node=int(i),psi_V=V*dx[i],phin_V=V*dx[N+i],phip_V=V*dx[2*N+i]))
        jv=[]
        for label in ('zero','finite_fixed'):jv+=c.s.prior.old.jc.jvp_metrics(root/label/'jvp.csv',case['key'],label)[0]
        edges0=a.rows(zero/'edges.csv');edges1=a.rows(changed/'edges.csv');unchanged=all(x[k]==y[k] for x,y in zip(edges0,edges1) for k in ('electron_flux','hole_flux','electron_mobility_m2_V_s','hole_mobility_m2_V_s','couple_m'))
        passed=np.array_equal(R,original) and source_error<=1e-5 and weak<=1e-12 and linear<=1e-8 and all(r['qualified'] for r in jv) and unchanged
        checks.append(dict(key=case['key'],zero_identity=bool(np.array_equal(R,original)),source_relative=source_error,weak_J_relative=weak,linear_relative=linear,transport_identity=unchanged,qualified=bool(passed)))
        print(checks[-1],flush=True)
    a.write_csv(OUT/'preflight.csv',checks);a.write_csv(OUT/'tangent.csv',tangents);d.matrix.freeze(OUT/'preflight_evidence.json',[OUT/'freeze.json',OUT/'preflight.csv',OUT/'tangent.csv']+[f for case in a.read(OUT/'contract.json')['cases'] for f in (LOCAL/case['key']/'fixed').rglob('*') if f.is_file()])


def run():
    a.verify(OUT/'preflight_evidence.json');assert all(x['qualified']=='True' for x in a.rows(OUT/'preflight.csv'))
    def one(case):
        rows=[]
        for job in case['jobs']:
            dest=Path(job['dest']);status=execute(dest/'config.json',case,job['alpha']);row=dict(key=case['key'],label=job['label'],qualified=False,failure=status.get('failure_reason',''))
            for name in ('all_row','acceptance_edges'):execute(dest/(name+'.json'),case,job['alpha'])
            row.update(c.s.prior.old.qualify(case,dest));rows.append(row);print(case['key'],job['label'],row['qualified'],row['max_row_ratio'],flush=True)
        return rows
    with ThreadPoolExecutor(max_workers=2) as pool:rows=[r for group in pool.map(one,a.read(OUT/'contract.json')['cases']) for r in group]
    keys=list(dict.fromkeys(k for r in rows for k in r));a.write_csv(OUT/'dc.csv',[{k:r.get(k,'') for k in keys} for r in rows])


def analyze():
    dc=a.rows(OUT/'dc.csv');dr=a.rows(OUT/'tangent.csv');out=[]
    for case in a.read(OUT/'contract.json')['cases']:
        geo,mask=v.m.previous.prior.support(case);ids=np.flatnonzero(mask);root=LOCAL/case['key'];base=d.ordered(Path(case['baseline'])/'state.csv',geo.count);zero=d.ordered(root/'zero/state.csv',geo.count)
        noise=float(np.linalg.norm(t.physical_delta(zero,base,ids)));delta={name:t.physical_delta(d.ordered(root/name/'state.csv',geo.count),zero,ids) for name in ('plus_full','minus_full','plus_half','minus_half')}
        odd={amp:(delta['plus_'+amp]-delta['minus_'+amp])/2 for amp in ('full','half')};tr=[r for r in dr if r['key']==case['key']];assert [int(r['node']) for r in tr]==list(ids)
        tangent=np.array([[float(r[f+'_V']) for r in tr] for f in ('psi','phin','phip')]);lin=float(np.linalg.norm(odd['full']-2*odd['half'])/np.linalg.norm(odd['full']))
        for amp,alpha in [('full',.001),('half',.0005)]:
            signal=float(np.linalg.norm(odd[amp]));error=float(np.linalg.norm(odd[amp]-alpha*tangent)/np.linalg.norm(alpha*tangent));even=float(np.linalg.norm((delta['plus_'+amp]+delta['minus_'+amp])/2)/signal);snr=signal/max(noise,1e-300)
            passed=all(r['qualified']=='True' for r in dc if r['key']==case['key'] and r['label'] in ('zero','plus_'+amp,'minus_'+amp)) and error<=1e-3 and lin<=1e-3 and even<=.01 and snr>=100
            out.append(dict(key=case['key'],amplitude=amp,prediction_relative=error,two_amplitude_relative=lin,even_over_odd=even,signal_over_zero_drift=snr,qualified=passed))
    a.write_csv(OUT/'response.csv',out);summary=dict(DC=len(dc),qualified_DC=sum(r['qualified']=='True' for r in dc),responses=len(out),qualified_responses=sum(r['qualified'] for r in out),qualified=all(r['qualified'] for r in out),native_response_qualified=False,finite_replacement=False)
    a.write(OUT/'summary.json',summary);d.matrix.freeze(OUT/'evidence.json',[OUT/'freeze.json',OUT/'preflight_evidence.json',OUT/'summary.json',OUT/'dc.csv',OUT/'response.csv']+[x for case in a.read(OUT/'contract.json')['cases'] for job in case['jobs'] for x in Path(job['dest']).rglob('*') if x.is_file()]);print(summary,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('build','prepare','preflight','run','analyze'));globals()[parser.parse_args().action]()
