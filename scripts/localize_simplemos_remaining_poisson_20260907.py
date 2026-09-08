"""Independent remaining-Poisson ledger at the two qualified eight-point states."""
import argparse, math
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
import complete_simplemos_production_migration_20260907 as v
import qualify_simplemos_native_mobility_cache_20260907 as q
import audit_simplemos_masetti_native_geometry_20260907 as g
import audit_simplemos_masetti_poisson_constants_20260907 as constants

a=v.a; d=v.d; p=v.p
LOCAL=p.REPO/'build-release/simplemos_remaining_poisson_20260907'
OUT=p.REPO/'reference_tcad/simplemos_sentaurus2022/remaining_poisson_20260907'

def env(c, mobility=1):
    e=v.environment()
    e.update(VELA_CANDIDATE_MOBILITY_ALPHA=str(mobility),
             VELA_CANDIDATE_MOBILITY_RATIOS=str(q.LOCAL/'ratios'/(c['device']+'.txt')))
    return e

def live(edges, count):
    values={}
    for e in edges:
        for side in ('0','1'):
            i=int(e['node'+side])
            value=tuple(float(e[k+side+u]) for k,u in [('psi','_V'),('electron_density','_m3'),('hole_density','_m3')])
            if i in values: assert values[i]==value
            values[i]=value
    assert len(values)==count
    return np.array([values[i] for i in range(count)]).T

def raw_root(c):
    base=Path(c['base'])
    return g.p.prev.LOCAL/'native_exports/masetti'/base.parent.name/f"vg_{int(base.name.split('_')[1]):03d}"

def prepare():
    a.verify(q.OUT/'final_evidence.json');a.verify(v.OUT/'validation_evidence.json')
    cases=a.read(v.OUT/'contract.json')['cases'];files=[Path(__file__).resolve(),q.OUT/'final_evidence.json',v.OUT/'validation_evidence.json',q.RUNNER]
    for c in cases:
        c['ledger_probes']=[]
        for label,mu,state in [('joint',0,v.LOCAL/c['key']/'from_vela/state.csv'),
                              ('native_mu',1,q.LOCAL/c['key']/'replacement/state.csv'),
                              ('native_state',1,v.m.previous.LOCAL/c['key']/'native_referenced_joint/state.csv')]:
            paths=v.probes(v.config(c),LOCAL/c['key']/label,state)
            c['ledger_probes'] += [dict(path=str(x),mobility=mu) for x in paths if x.stem in ('functional','edges') or (x.stem=='adjoint' and label!='native_state')]
            files+=paths+[state]
    a.write(OUT/'ledger_contract.json',dict(cases=cases,scope='Read-only live-state Poisson residual and signed adjoint screening; no response calibration inferred from earlier different directions.',
        components=['native replay','native geometry','epsilon0 convention','charge convention','density mapping','packing','remainder'],acceptance_changed=False))
    d.matrix.freeze(OUT/'ledger_freeze.json',files+[OUT/'ledger_contract.json'])

def run():
    a.verify(OUT/'ledger_freeze.json')
    def one(c):
        for j in c['ledger_probes']:
            s=v.execute(Path(j['path']),q.RUNNER,env(c,j['mobility']))
            assert s['exit_code']==0,(j,s)
        print('Poisson ledger probes',c['key'],flush=True)
    with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(one,a.read(OUT/'ledger_contract.json')['cases']))

def analyze():
    a.verify(OUT/'ledger_freeze.json')
    factor=a.read(d.matrix.p.CONTRACT)['frozen_poisson_residual_per_physical_charge']
    cases=a.read(OUT/'ledger_contract.json')['cases'];checks=[];ledger=[];spatial=[];top=[];geometry=[]
    for device in ('n19','n23'):
        geo,xy,elements,volumes,K,parts,info,_=g.geometry(device)
        nativeK=K*(constants.NATIVE_EPS0/d.matrix.spatial.m73.EPS0);legacy=geo.matrices['legacy']
        geometry.append(info)
        for c in (c for c in cases if c['device']==device):
            root=LOCAL/c['key'];psi,n,h,spread=d.matrix.spatial.m73.sentaurus_state(raw_root(c),geo);assert spread<=1e-12
            rows=d.ordered(root/'native_state/residual.csv',geo.count)
            net=np.array([float(x['net_doping_m3']) for x in rows])
            charge=n-h-net;packed=live(a.rows(root/'native_state/edges.csv'),geo.count)
            nativeR=nativeK@psi+constants.NATIVE_Q*charge*volumes['Si']
            sources={
                'native_constant_replay':nativeR,
                'dielectric_geometry':(legacy-K)@psi,
                'epsilon0_convention':(K-nativeK)@psi,
                'charge_constant_convention':(d.fixed.Q-constants.NATIVE_Q)*charge*volumes['Si'],
                'thermal_density_mapping':d.fixed.Q*((packed[1]-n)-(packed[2]-h))*volumes['Si'],
                'potential_packing':legacy@(packed[0]-psi)}
            actual=np.array([float(x['psi_residual']) for x in rows])/factor
            sources['replay_remainder']=actual-sum(sources.values())
            den=np.linalg.norm(actual[geo.free]);charge_norm=np.linalg.norm((constants.NATIVE_Q*charge*volumes['Si'])[geo.free])
            check=dict(key=c['key'],native_replay_over_charge=float(np.linalg.norm(nativeR[geo.free])/charge_norm),
                       live_replay_relative=float(np.linalg.norm(sources['replay_remainder'][geo.free])/den))
            check['qualified']=check['native_replay_over_charge']<=1e-10 and check['live_replay_relative']<=1e-8
            checks.append(check)
            for label in ('joint','native_mu'):
                weights=d.ordered(root/label/'adjoint.csv',geo.count)
                lam=np.array([float(x['lambda_poisson']) for x in weights]);lam[geo.contact_nodes]=0
                Id=a.read(root/label/'functional.status.json')['current_A_per_um']
                def project(source):return -math.fsum(float(x*y) for x,y in zip(lam,source*factor))
                for name,source in sources.items():
                    value=project(source)
                    ledger.append(dict(key=c['key'],state=label,component=name,mapped_residual_projection_A_per_um=value,relative_to_Id=value/Id,
                                       residual_norm_C_per_m=float(np.linalg.norm(source[geo.free])),screening_only=True))
                # Candidate derivative is evaluated at the qualified state, not the native mapped state.
                state=live(a.rows(root/label/'edges.csv'),geo.count)
                groups={};edges=[]
                for (i,j),part in parts.items():
                    mats=';'.join(sorted({x['material'] for x in part}))
                    delta=float(-K[i,j]+legacy[i,j]);s=delta*(state[0,i]-state[0,j])
                    value=-(lam[i]-lam[j])*s*factor
                    groups.setdefault(mats,[]).append(value)
                    edges.append(dict(key=c['key'],state=label,node0=i,node1=j,materials=mats,delta_K_F_per_m=delta,
                                      conservative_source_C_per_m=s,current_derivative_A_per_um=value,relative_to_Id=value/Id))
                for name,values in groups.items():
                    spatial.append(dict(key=c['key'],state=label,materials=name,edges=len(values),signed_prediction_A_per_um=math.fsum(values),
                                        absolute_contributions_A_per_um=math.fsum(abs(x) for x in values),relative_to_Id=math.fsum(values)/Id,screening_only=True))
                top+=sorted(edges,key=lambda x:abs(x['current_derivative_A_per_um']),reverse=True)[:30]
                # Native-convention changes are separate axes, never fitted to Id.
                for name,source in [('candidate_dielectric_geometry',(K-legacy)@state[0]),
                                    ('candidate_epsilon0',(constants.NATIVE_EPS0/d.matrix.spatial.m73.EPS0-1)*(legacy@state[0])),
                                    ('candidate_charge_constant',(constants.NATIVE_Q-d.fixed.Q)*(state[1]-state[2]-net)*volumes['Si'])]:
                    value=project(source)
                    ledger.append(dict(key=c['key'],state=label,component=name,mapped_residual_projection_A_per_um=value,relative_to_Id=value/Id,
                                       residual_norm_C_per_m=float(np.linalg.norm(source[geo.free])),screening_only=True))
    for name,rows in [('ledger_checks',checks),('ledger',ledger),('dielectric_spatial',spatial),('dielectric_top_edges',top),('geometry',geometry)]:a.write_csv(OUT/(name+'.csv'),rows)
    d.matrix.freeze(OUT/'ledger_evidence.json',[OUT/'ledger_freeze.json']+list(OUT.glob('*.csv'))+[x for c in cases for x in (LOCAL/c['key']).rglob('*') if x.is_file()])
    print('Live and native Poisson replays qualified:',sum(x['qualified'] for x in checks),'/',len(checks),flush=True)
    print([x for x in ledger if x['state']=='native_mu' and x['component'].startswith('candidate_')],flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('prepare','run','analyze'));globals()[parser.parse_args().action]()
