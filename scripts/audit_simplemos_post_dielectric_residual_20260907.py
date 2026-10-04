"""Actual-state first-order ledger after K qualification, including boundary rows.

This is not a same-source calibration of a new physical-constant candidate.
The mapped native state may violate Vela boundary conventions. Its boundary
rows must be included when comparing states, unlike a fixed-boundary K source.
"""
import argparse, math
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import validate_simplemos_remaining_dielectric_20260907 as t

a=t.a;d=t.d;l=t.l;OUT=t.OUT/'post_K_ledger';LOCAL=t.LOCAL/'post_K_ledger'

def prepare():
    a.verify(t.OUT/'validation_evidence.json');assert a.read(t.OUT/'summary.json')['all_native_K_promoted']
    cases=a.read(t.OUT/'contract.json')['cases'];files=[Path(__file__).resolve(),t.OUT/'validation_evidence.json',t.RUNNER]
    for c in cases:
        c['post_K_probes']=[]
        for label,state in [('qualified',t.LOCAL/c['key']/'all_native_K/replacement/state.csv'),('native',t.v.m.previous.LOCAL/c['key']/'native_referenced_joint/state.csv')]:
            paths=t.v.probes(t.v.config(c),LOCAL/c['key']/label,state);files+=paths+[state]
            c['post_K_probes'] += [str(x) for x in paths if (label=='qualified' and x.stem=='adjoint') or (label=='native' and x.stem in ('functional','edges','terms'))]
    a.write(OUT/'contract.json',dict(cases=cases,read_only=True,interpretation='First-order actual-state residual ledger, including contact boundary rows. New q/epsilon0/SRH source directions remain uncalibrated.',
        identity='Id_strict-Id_native = -lambda dot (R_mapped-R_strict) + (Id_mapped-Id_native) + finite_state_remainder',
        fitted_parameters=False,acceptance_changed=False))
    d.matrix.freeze(OUT/'freeze.json',files+[OUT/'contract.json'])

def run():
    a.verify(OUT/'freeze.json')
    def one(c):
        for path in c['post_K_probes']:
            s=t.execute(Path(path),c,'all_native_K',1.);assert s['exit_code']==0,(path,s)
        print('Post-K residual probes',c['key'],flush=True)
    with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(one,a.read(OUT/'contract.json')['cases']))

def analyze():
    a.verify(OUT/'freeze.json');ledger=[];checks=[];factor=a.read(d.matrix.p.CONTRACT)['frozen_poisson_residual_per_physical_charge']
    for c in a.read(OUT/'contract.json')['cases']:
        root=LOCAL/c['key'];strict=t.LOCAL/c['key']/'all_native_K/replacement/post'
        geo,xy,el,vol,K,parts,*_=l.g.geometry(c['device']);free=geo.free;contacts=geo.contact_nodes
        nativeK=K*l.constants.NATIVE_EPS0/d.matrix.spatial.m73.EPS0
        Rn=d.array(d.ordered(root/'native/residual.csv',geo.count),('psi_residual','phin_residual','phip_residual'))
        Rs=d.array(d.ordered(strict/'residual.csv',geo.count),('psi_residual','phin_residual','phip_residual'))
        lam=d.array(d.ordered(root/'qualified/adjoint.csv',geo.count),('lambda_poisson','lambda_electron','lambda_hole'))
        psi,n,h,spread=d.matrix.spatial.m73.sentaurus_state(l.raw_root(c),geo);assert spread<=1e-12
        packed=l.live(a.rows(root/'native/edges.csv'),geo.count)
        net=np.array([float(r['net_doping_m3']) for r in d.ordered(root/'native/residual.csv',geo.count)])
        charge=n-h-net
        poisson={'Poisson_native_replay':nativeK@psi+l.constants.NATIVE_Q*charge*vol['Si'],
            'Poisson_epsilon0_convention':(K-nativeK)@psi,
            'Poisson_charge_constant_convention':(d.fixed.Q-l.constants.NATIVE_Q)*charge*vol['Si'],
            'Poisson_thermal_density_mapping':d.fixed.Q*((packed[1]-n)-(packed[2]-h))*vol['Si'],
            'Poisson_potential_packing':K@(packed[0]-psi)}
        poisson['Poisson_replay_remainder']=Rn[0]/factor-sum(poisson.values())
        poisson['Poisson_minus_strict_residual']=-Rs[0]/factor
        components={}
        for name,source in poisson.items():
            value=np.zeros_like(Rn);value[0,free]=source[free]*factor;components[name]=value
        nt=d.ordered(root/'native/terms.csv',geo.count);st=d.ordered(strict/'terms.csv',geo.count)
        for block,car in [(1,'electron'),(2,'hole')]:
            for field in ('flux','recombination'):
                value=np.zeros_like(Rn);value[block,free]=np.array([float(nt[i][car+'_'+field])-float(st[i][car+'_'+field]) for i in free])
                components[car+'_'+('SRH_state_difference' if field=='recombination' else 'flux_state_difference')]=value
        boundary=np.zeros_like(Rn);boundary[:,contacts]=(Rn-Rs)[:,contacts];components['contact_boundary_state_difference']=boundary
        components['continuity_other_and_roundoff']=Rn-Rs-sum(components.values())
        projected=lambda source:-math.fsum(float(x*y) for x,y in zip(lam.flat,source.flat))
        actual=t.a.read(t.LOCAL/c['key']/'all_native_K/replacement/config.status.json')['contact_currents_A_per_um']['drain']-c['native_Id_A_per_um']
        direct=a.read(root/'native/functional.status.json')['current_A_per_um']-c['native_Id_A_per_um']
        pred=math.fsum(projected(x) for x in components.values())+direct
        for name,source in components.items():
            value=projected(source)
            ledger.append(dict(key=c['key'],component=name,current_A_per_um=value,relative_to_native_Id=value/c['native_Id_A_per_um'],screening_only=True))
        for name,value in [('mapped_port_minus_native_port',direct),('finite_state_remainder',actual-pred),('actual_remaining_Id_error',actual)]:
            ledger.append(dict(key=c['key'],component=name,current_A_per_um=value,relative_to_native_Id=value/c['native_Id_A_per_um'],screening_only=name!='actual_remaining_Id_error'))
        closure=float(np.linalg.norm(Rn-Rs-sum(components.values()))/max(np.linalg.norm(Rn-Rs),1e-300))
        checks.append(dict(key=c['key'],source_accounting_relative=closure,actual_remaining_A_per_um=actual,first_order_sum_A_per_um=pred,
            finite_state_remainder_over_Id=(actual-pred)/c['native_Id_A_per_um'],finite_state_remainder_over_gap=(actual-pred)/actual,
            new_parameter_response_calibrated=False))
    a.write_csv(OUT/'ledger.csv',ledger);a.write_csv(OUT/'closure.csv',checks)
    d.matrix.freeze(OUT/'evidence.json',[OUT/'freeze.json',OUT/'ledger.csv',OUT/'closure.csv']+[x for x in LOCAL.rglob('*') if x.is_file()])
    print('Post-K first-order ledger:',checks,flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('prepare','run','analyze'));globals()[parser.parse_args().action]()
