"""Ascending strict Vela continuation with explicit failure stop and per-state audit."""
import argparse
import copy
from concurrent.futures import ThreadPoolExecutor
import math
from pathlib import Path
import validate_simplemos_local_conservative_flux as v
import prepare_simplemos_fullfield_native_20260906 as native

a=v.a;d=v.prior;LOCAL=native.LOCAL;OUT=native.OUT


def prepare():
    a.verify(v.FREEZE);a.verify(v.OUT/'evidence.json')
    files=[Path(__file__).resolve(),v.RUNNER,v.FREEZE,v.OUT/'evidence.json']
    cases=[]
    for c in a.read(d.CONTRACT)['cases']:
        cfg=a.read(d.matrix.spatial.precision.LOCAL/c['case']/'config.json')
        initial=Path(cfg['state_file']).parent.parent/'10_drain_ramp/state.csv'
        assert initial.exists(),initial
        files.append(initial)
        cfg.pop('simplemos_m3',None)
        paths=[]
        for index in range(21):
            dest=LOCAL/'vela'/c['case']/f'vg_{index:03d}'
            point=copy.deepcopy(cfg)
            for contact in point['contacts']:
                if contact['name']=='gate':contact['bias']=round(index*.05,12)
            point['state_file']=str(initial if index==0 else dest.parent/f'vg_{index-1:03d}'/'state.csv')
            point['output_state_file']=str(dest/'state.csv')
            path=dest/'config.json';a.write(path,point);files.append(path);paths.append(str(path))
        files += [Path(cfg[k]) for k in ('mesh_file','materials_file','node_doping_file')]
        cases.append(dict(case=c['case'],device=c['device'],vd=float(c['vd']),configs=paths))
    a.write(OUT/'vela_contract.json',dict(status='frozen_before_execution',cases=cases,
        solver='Frozen isolated vector-chain corrected runner; original physical residual, SparseLU, strict carrier acceptance unchanged.',
        continuation='Ascending Vg .05 V increments from existing M65 drain state; stop a branch after first failed strict acceptance. No continuation from a failed state.',
        recovery='If needed, separately freeze native-coherent initial guesses for failed/unreached target points; keep ascending failures; unchanged physical model and strict acceptance.',
        gates={'carrier_eps_row':1e-6,'global_tolerance':1e-6,'global_source_floor':1e-10,'kcl_over_Id':1e-8},
        possible_target_solves=84,m82_released=False,m83_released=False))
    files.append(OUT/'vela_contract.json');d.matrix.freeze(OUT/'vela_freeze.json',files)
    print('Frozen 84 Vela ascending target configurations',flush=True)


def qualify(dest):
    cfgpath=dest/'config.json';s=v.execute(cfgpath)
    row=dict(case=dest.parent.name,vg=int(dest.name[-3:])*.05,exit_code=s['exit_code'],converged=s.get('converged',False),
        reason=s.get('convergence_reason',''),iterations=s.get('iterations',0),qualified=False,current_A_per_um=s.get('contact_currents_A_per_um',{}).get('drain',0.))
    if (dest/'state.csv').exists():
        cfg=a.read(cfgpath);cfg.pop('output_state_file');cfg.update(state_file=str(dest/'state.csv'),simulation_type='newton_carrier_term_probe',output_csv=str(dest/'carrier.csv'),carrier_term_probe={'solved_equation_terms':True})
        cfg['solver']['carrier_row_convergence']['mode']='report';cfg['solver']['global_continuity_closure']={'mode':'enforce','tolerance':1e-6,'source_floor':1e-10}
        if not (dest/'acceptance.json').exists():a.write(dest/'acceptance.json',cfg)
        audit=v.execute(dest/'acceptance.json')
        cc=s.get('contact_currents_A_per_um',{})
        kcl=abs(math.fsum(cc.values()))/max(abs(row['current_A_per_um']),1e-300)
        row.update(kcl_over_Id=kcl,local_violations=audit['carrier_row_convergence']['violation_count'],
            global_electron_qualified=audit['global_continuity_closure']['electron']['qualified'],global_hole_qualified=audit['global_continuity_closure']['hole']['qualified'])
        row['qualified']=s['exit_code']==audit['exit_code']==0 and s['converged'] and audit['carrier_row_convergence']['satisfied'] and audit['global_continuity_closure']['satisfied'] and kcl<=1e-8
    return row


def curve(case):
    ledger=[]
    for path in case['configs']:
        dest=Path(path).parent
        inputpath=Path(a.read(Path(path))['state_file'])
        if not (dest/'state_input_freeze.json').exists():d.matrix.freeze(dest/'state_input_freeze.json',[inputpath,Path(path)])
        a.verify(dest/'state_input_freeze.json')
        row=qualify(dest);ledger.append(row)
        print(row,flush=True)
        if not row['qualified']:break
    if not (LOCAL/'vela'/case['case']/'ledger.json').exists():a.write(LOCAL/'vela'/case['case']/'ledger.json',ledger)
    return ledger


def run():
    a.verify(OUT/'vela_freeze.json')
    with ThreadPoolExecutor(max_workers=2) as pool:groups=list(pool.map(curve,a.read(OUT/'vela_contract.json')['cases']))
    a.write_csv(OUT/'vela_ascending.csv',[r for g in groups for r in g])


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('prepare','run'))
    globals()[parser.parse_args().action]()
