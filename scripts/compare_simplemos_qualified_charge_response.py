"""Screen spatial charge responses on strictly qualified states and M81 weights."""
import argparse
import json
import math
import numpy as np
from pathlib import Path
import run_simplemos_strict_flux_precision as precision
import run_simplemos_fixed_state_formula_audit as fixed
import audit_simplemos_m81_native_green_integral as native

old=fixed.old;m73=fixed.m73
REPO,ROOT=fixed.REPO,fixed.ROOT
LOCAL=REPO/'build-release/simplemos_qualified_charge_response'
OUT=ROOT/'qualified_charge_response'
CONTRACT=OUT/'contract.json'
SCRIPT=Path(__file__).resolve()


def native_arrays():
    cfg=old.read_json(native.CONTRACT)
    measures=native.m34.parse_debug_block(native.DEBUG.read_text(),'Measure')
    green={int(r['node_id']):float(r['component0']) for r in native.csvrows(native.LOCAL/'green_exports/n23/fields/CurPotReACGreenFunction_region0.csv')}
    v={i:0. for i in green}
    for e in native.csvrows(native.M34/'n23_import/elements.csv'):
        if e['material']!='Si':continue
        for j,k in enumerate(cfg['measure_permutation']):v[int(e[f'node{j}'])]+=measures[int(e['id'])]['values'][k]
    return v,green


def prepare():
    if CONTRACT.exists():raise FileExistsError(CONTRACT)
    precision.audit.verify(precision.FREEZE)
    cases=old.read_csv(precision.OUT/'case_ledger.csv')
    if len(cases)!=4 or any(c['qualified']!='True' for c in cases):raise ValueError('Four strict states required')
    inputs=[SCRIPT,precision.FREEZE,precision.OUT/'case_ledger.csv',native.CONTRACT,native.FREEZE,native.DEBUG,
            native.M34/'n23_import/elements.csv',native.LOCAL/'green_exports/n23/fields/CurPotReACGreenFunction_region0.csv',
            native.OUT/'m81_native_pilot_result.json',fixed.upstream.original.upstream.previous.RUNNER]
    for c in cases:
        src=precision.LOCAL/c['case'];dest=LOCAL/c['case']
        deck=old.read_json(src/'config.json');deck.pop('output_state_file',None)
        deck.update(state_file=str(src/'state.csv'),simulation_type='terminal_current_adjoint_probe',contact='drain',
                    output_csv=str(dest/'adjoint.csv'),electron_volume_response={'output_csv':str(dest/'response.csv')})
        old.write_json(dest/'config.json',deck)
        inputs += [src/'config.json',src/'state.csv',dest/'config.json']
    old.write_json(CONTRACT,{'status':'frozen_before_execution','cases':cases,
        'source':'Positive fixed charge number density 1e13 cm^-3 added only to Poisson; unchanged mobility, dopants, SRH. source_F=-factor*q*1e19*area_m2. Response=-lambda dot source_F.',
        'supports':['uniform_Si','channel_0p05um','gate_interface_nodes'],
        'support_semantics':'x is depth; channel: Si within 0.05 um of surface and abs(y)<=0.125 um. Interface case selects Si/SiO2 shared nodes under gate; uses bulk nodal control volumes, not sheet density. Exclude Dirichlet Poisson nodes.',
        'geometry':'n23 uses authoritative M34 native signed Si areas for both codes; n19 uses geometric signed Si areas in Vela only. Uniform n23 native integral must recover M81 dI within existing 5% gate before spatial screening.',
        'scope':'Full-DD linear response, not a new self-consistent correction. M81 native comparison exists only for n23 low Vd. Localized native directions have no independent finite-difference qualification; they remain screening results even if the uniform calibration passes.',
        'gates':{'adjoint_relative_residual':1e-10,'current_identity_relative':1e-8,'poisson_scale_relative':1e-6,'uniform_native_integral_relative':.05},
        'new_nonlinear_solves':0,'new_sentaurus_runs':0,'m82_released':False,
        'inputs':{old.portable(p):old.sha256(p) for p in sorted(set(inputs))}})


def verify():
    for p,h in old.read_json(CONTRACT)['inputs'].items():
        if old.sha256(REPO/p)!=h:raise ValueError(f'Changed input {p}')


def run():
    verify();cfg=old.read_json(CONTRACT);rows=[];checks=[]
    old.m74.RUNNER=fixed.upstream.original.upstream.previous.RUNNER
    nv,green=native_arrays()
    for c in cfg['cases']:
        dest=LOCAL/c['case'];deck=old.read_json(dest/'config.json')
        status=old.run_probe(deck,dest/'config.json')
        geo=m73.Geometry(c['device']);masks,_,xy=old.m78.supports(c['device'],geo,.05)
        states={int(r['node_id']):r for r in old.read_csv(Path(deck['state_file']))}
        n=np.array([float(states[i]['electrons_m3']) for i in range(geo.count)])
        response=old.read_csv(dest/'response.csv');known=np.array([float(r['source_poisson']) for r in response])
        physical=m73.Q*n*(geo.volumes['barycentric_si']-geo.volumes['all_cell'])
        active=geo.free[np.abs(physical[geo.free])>np.max(np.abs(physical[geo.free]))*1e-6]
        factor=float(np.median(known[active]/physical[active]));factorerr=float(np.max(np.abs(known[active]/physical[active]/factor-1)))
        identity=abs(status['current_A_per_um']/float(c['current_A_per_um'])-1)
        if identity>1e-8 or factorerr>1e-6 or status['adjoint_relative_residual']>1e-10:raise ValueError('Response numerical gate failed')
        adj={int(r['node_id']):r for r in old.read_csv(dest/'adjoint.csv')}
        lam=np.array([float(adj[i]['lambda_poisson']) for i in range(geo.count)])
        volume=geo.volumes['signed_si'].copy()
        if c['device']=='n23':volume=np.array([nv.get(i,0)*1e-12 for i in range(geo.count)])
        free=np.zeros(geo.count,dtype=bool);free[geo.free]=True
        for name,mask in (('uniform_Si',masks['all_si']),('channel_0p05um',masks['channel']),('gate_interface_nodes',masks['gate_interface'])):
            mask=mask&free
            source=-factor*m73.Q*1e19*volume*mask
            vela=-float(lam@source)
            sent=None
            if c['device']=='n23' and float(c['vd'])==.05:
                sent=m73.Q*1e13*1e-12*math.fsum(green[i]*nv[i] for i in green if mask[i])
                if name=='uniform_Si':
                    reference=next(r['native_dI_A'] for r in old.read_json(native.OUT/'m81_native_pilot_result.json')['cases'] if r['device']=='n23')
                    err=abs(sent/reference-1)
                    if err>.05:raise ValueError('Uniform native control failed')
                    checks.append({'check':'n23_uniform_native_control','relative_error':err})
            rows.append({'case':c['case'],'device':c['device'],'vd':c['vd'],'support':name,'nodes':int(mask.sum()),
                'vela_delta_Id_A_per_um':vela,'vela_delta_log10_Id_linear':vela/(float(c['current_A_per_um'])*math.log(10)),
                'sentaurus_green_delta_Id_A_per_um':sent,'vela_to_native_response_ratio':None if sent is None else vela/sent,
                'native_direction_fd_qualified':name=='uniform_Si' and sent is not None})
        checks.append({'check':c['case']+'_adjoint_relative_residual','relative_error':status['adjoint_relative_residual']})
        print(c['case'],'response completed',flush=True)
    old.write_csv(OUT/'response_ledger.csv',rows);old.write_csv(OUT/'checks.csv',checks)
    old.write_json(OUT/'summary.json',{'strict_states':4,'vela_spatial_responses':12,'native_screened_directions':3,
        'native_independent_FD_qualified_directions':1,'m82_released':False,'new_nonlinear_solves':0,'new_sentaurus_runs':0})


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('prepare','run','verify'))
    {'prepare':prepare,'run':run,'verify':verify}[p.parse_args().action]()
