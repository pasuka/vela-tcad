"""Poisson reaction ledger along independently calibrated channel continuity response."""
import argparse
from decimal import Decimal
from pathlib import Path
import numpy as np
import validate_simplemos_local_conservative_flux as v

a=v.a;d=v.prior
LOCAL=d.REPO/'build-release/simplemos_channel_poisson_coupling_20260906'
OUT=d.ROOT/'channel_poisson_coupling_20260906'


def prepare():
    a.verify(v.FREEZE);a.verify(v.OUT/'evidence.json')
    files=[Path(__file__).resolve(),v.FREEZE,v.OUT/'evidence.json',v.RUNNER]
    cases=[]
    for c in a.read(v.CONTRACT)['cases']:
        r=next(r for r in c['regions'] if r['name']=='channel')
        for label,mult in (('plus_full',1.),('minus_full',-1.),('plus_half',.5),('minus_half',-.5)):
            source=v.LOCAL/c['case']/'channel'/label
            cfg=a.read(source/'config.json');cfg.pop('output_state_file')
            dest=LOCAL/c['case']/label
            cfg.update(simulation_type='newton_residual_probe',state_file=str(source/'state.csv'),output_csv=str(dest/'residual.csv'))
            path=dest/'config.json';a.write(path,cfg)
            files += [path,source/'state.csv',source/'config.status.json',source/'acceptance.status.json']
            cases.append(dict(case=c['case'],label=label,edge=r['edge'],native_flux=mult*r['native_amplitude'],config=str(path)))
    a.write(OUT/'contract.json',dict(status='frozen_before_execution',cases=cases,
        definition='At equal bias, decompose delta psi into K^-1 Poisson residual change, electron charge, hole charge and boundary response; use actual +/- strict channel-current states. K and charge volumes are the unchanged Vela operator.',
        gates={'poisson_response_relative':1e-3,'full_half_field_relative':1e-3},
        limitations='Conditional charge-response bookkeeping, not an independent Sentaurus Poisson operator or a qualified physical correction.',new_nonlinear_solves=0,new_native_runs=0))
    files.append(OUT/'contract.json');d.matrix.freeze(OUT/'freeze.json',files)


def run():
    a.verify(OUT/'freeze.json')
    for c in a.read(OUT/'contract.json')['cases']:
        s=v.execute(Path(c['config']),c['edge'],c['native_flux']);assert s['exit_code']==0
    print('Eight calibrated-state Poisson residual probes complete',flush=True)


def analyze():
    a.verify(OUT/'freeze.json');ledger=[];checks=[];nodeledger=[]
    for c in a.read(v.CONTRACT)['cases']:
        key=c['case'];count=c['nodes'];geo=d.matrix.spatial.m73.Geometry('n23')
        masks,_,xy=d.matrix.spatial.old.m78.supports('n23',geo,.05)
        K=geo.matrices['legacy'];volume=geo.volumes['all_cell']
        jobs=[]
        for amplitude in ('full','half'):
            states=[d.ordered(v.LOCAL/key/'channel'/(sign+'_'+amplitude)/'state.csv',count) for sign in ('plus','minus')]
            delta=np.array([[float((Decimal(p[col])-Decimal(m[col]))/2) for p,m in zip(*states)] for col in ('psi','electrons_m3','holes_m3')])
            residuals=[d.ordered(LOCAL/key/(sign+'_'+amplitude)/'residual.csv',count) for sign in ('plus','minus')]
            dr=np.array([float(p['psi_residual'])-float(m['psi_residual']) for p,m in zip(*residuals)])/2/c['factor']
            jobs.append((amplitude,delta,dr))
        states=[d.ordered(d.matrix.spatial.precision.LOCAL/key/'state.csv',count),d.ordered(Path(c['mapped']),count)]
        delta=d.array(states[0],('psi','electrons_m3','holes_m3'))-d.array(states[1],('psi','electrons_m3','holes_m3'))
        rr=[d.ordered(d.LOCAL/key/role/'residual.csv',count) for role in ('strict','mapped')]
        dr=np.array([float(p['psi_residual'])-float(m['psi_residual']) for p,m in zip(*rr)])/c['factor']
        jobs.append(('absolute_strict_minus_mapped',delta,dr))
        responses={}
        for name,delta,dr in jobs:
            components={}
            for carrier,source in (('electron',d.fixed.Q*delta[1]*volume),('hole',-d.fixed.Q*delta[2]*volume)):
                for region in ('channel','source','drain','substrate'):
                    components[carrier+'_'+region]=geo.solve('legacy',source*masks[region])
            boundary=np.zeros(count);boundary[geo.contact_nodes]=delta[0,geo.contact_nodes]
            boundary+=geo.solve('legacy',K@boundary)
            components['boundary']=boundary
            components['operator_residual']=geo.solve('legacy',-dr)
            predicted=sum(components.values())
            norm=np.linalg.norm(delta[0][geo.free]);error=np.linalg.norm((predicted-delta[0])[geo.free])/max(norm,1e-300)
            charge_only=sum(value for k,value in components.items() if k.startswith(('electron','hole')))+boundary
            closure=np.linalg.norm((charge_only-delta[0])[geo.free])/max(norm,1e-300)
            checks.append(dict(case=key,comparison=name,full_reconstruction_relative=float(error),charge_only_relative=float(closure),
                charge_response_qualified=name!='absolute_strict_minus_mapped' and closure<=1e-3,
                maximum_actual_psi_V=float(max(abs(delta[0])))))
            responses[name]=delta[0]
            for component,value in components.items():
                ledger.append(dict(case=key,comparison=name,component=component,
                    at_interface_338_V=float(value[338]),at_channel_320_V=float(value[320]),at_channel_324_V=float(value[324]),
                    channel_rms_V=float(np.sqrt(np.average(value[masks['channel']]**2,weights=geo.volumes['barycentric_si'][masks['channel']])))))
            for node in (320,324,338):
                nodeledger.append(dict(case=key,comparison=name,node_id=node,actual_psi_V=float(delta[0,node]),reconstructed_psi_V=float(predicted[node])))
        linearity=float(np.linalg.norm(responses['full']-2*responses['half'])/np.linalg.norm(responses['full']))
        for row in checks:
            if row['case']==key:row['full_half_psi_relative']=linearity
    for name,rows in (('components',ledger),('checks',checks),('nodes',nodeledger)):a.write_csv(OUT/(name+'.csv'),rows)
    print(checks,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('prepare','run','analyze'))
    globals()[parser.parse_args().action]()
