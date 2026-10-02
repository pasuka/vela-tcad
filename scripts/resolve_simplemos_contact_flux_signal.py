"""New, separately frozen signal-resolution control; preserve original failed gates."""
import argparse
import copy
from decimal import Decimal
import math
import numpy as np
import validate_simplemos_local_conservative_flux as v

a=v.a
LOCAL=v.LOCAL/'contact_resolution'
OUT=v.OUT/'contact_resolution'
FREEZE=OUT/'freeze.json'
CONTRACT=OUT/'contract.json'


def prepare():
    a.verify(v.FREEZE)
    cases=[];files=[__import__('pathlib').Path(__file__).resolve(),v.FREEZE,v.RUNNER]
    for c0 in a.read(v.CONTRACT)['cases']:
        c=copy.deepcopy(c0)
        r=next(r for r in c['regions'] if r['name']=='drain_contact')
        gain=r['total_A_per_um']/r['current_amplitude_A_per_um']
        relative=10.**math.ceil(math.log10(1e-7/abs(gain)))
        c.update(relative_amplitude=relative,multiplier=relative/r['relative_current'])
        deck=a.read(v.LOCAL/c['case']/'zero/config.json')
        for label,mult in (('plus_full',1.),('minus_full',-1.),('plus_half',.5),('minus_half',-.5)):
            dest=LOCAL/c['case']/label
            d=copy.deepcopy(deck);d['output_state_file']=str(dest/'state.csv')
            a.write(dest/'config.json',d);files.append(dest/'config.json')
        cases.append(c)
    a.write(CONTRACT,{'status':'frozen_before_execution','cases':cases,
        'amplitude_rule':'Smallest power of ten in units of baseline Id whose frozen predicted net response reaches 1e-7 Id. Both signs and half amplitude; no fit to new DC values.',
        'gates':a.read(v.CONTRACT)['gates'],'additional_max_state_change_V':.001,
        'original_small_amplitude_evidence_retained':True,'new_nonlinear_solves':8,'new_sentaurus_runs':0})
    files.append(CONTRACT);v.m.freeze(FREEZE,files)
    print('Frozen separate contact resolution control',[(c['vd'],c['relative_amplitude']) for c in cases],flush=True)


def run():
    a.verify(FREEZE);ledger=[]
    for c in a.read(CONTRACT)['cases']:
        region=next(r for r in c['regions'] if r['name']=='drain_contact')
        for label,mult in (('plus_full',1.),('minus_full',-1.),('plus_half',.5),('minus_half',-.5)):
            dest=LOCAL/c['case']/label
            native=region['native_amplitude']*c['multiplier']*mult
            s=v.execute(dest/'config.json',region['edge'],native)
            d=a.read(dest/'config.json');d.pop('output_state_file');d['state_file']=str(dest/'state.csv')
            d['solver']['carrier_row_convergence']['mode']='report'
            d['solver']['global_continuity_closure']=a.read(v.CONTRACT)['global_profile']
            d.update(simulation_type='newton_carrier_term_probe',output_csv=str(dest/'terms.csv'),carrier_term_probe={'solved_equation_terms':True})
            a.write(dest/'acceptance.json',d)
            audit=v.execute(dest/'acceptance.json',region['edge'],native)
            cc=s['contact_currents_A_per_um'];cur=cc['drain'];kcl=abs(math.fsum(cc.values()))/abs(cur)
            source=a.read(dest/'config.json')['state_file']
            initial=v.prior.ordered(__import__('pathlib').Path(source),c['nodes'])
            final=v.prior.ordered(dest/'state.csv',c['nodes'])
            changes=[]
            for u,w in zip(initial,final,strict=True):
                changes.append(abs(float(w['psi'])-float(u['psi'])))
                for carrier in ('electron','hole'):
                    ref,inc=carrier+'_qf_reference_V',carrier+'_qf_increment_V'
                    changes.append(abs(float(Decimal(w[ref])-Decimal(u[ref])+Decimal(w[inc])-Decimal(u[inc]))))
            maximum=max(changes)
            row={'case':c['case'],'label':label,'current_A_per_um':cur,'iterations':s['iterations'],'reason':s['convergence_reason'],
                'max_state_change_V':maximum,'kcl_over_Id':kcl,'local_violations':audit['carrier_row_convergence']['violation_count'],
                'qualified':s['exit_code']==0 and audit['exit_code']==0 and s['converged'] and audit['carrier_row_convergence']['satisfied'] and audit['global_continuity_closure']['satisfied'] and kcl<=1e-8 and maximum<=.001}
            ledger.append(row);print(row,flush=True)
    a.write_csv(OUT/'dc.csv',ledger)


def analyze():
    a.verify(FREEZE);rows=a.rows(OUT/'dc.csv');cal=[]
    for c in a.read(CONTRACT)['cases']:
        region=next(r for r in c['regions'] if r['name']=='drain_contact')
        group=[r for r in rows if r['case']==c['case']]
        values={r['label']:float(r['current_A_per_um']) for r in group}
        zero=a.read(v.LOCAL/c['case']/'zero/config.status.json')['contact_currents_A_per_um']['drain']
        drift=abs(zero-float(c['current_A_per_um']))
        full=(values['plus_full']-values['minus_full'])/2;half=(values['plus_half']-values['minus_half'])/2
        linearity=abs(2*half/full-1)
        for label,mult,fd in (('full',1.,full),('half',.5,half)):
            predicted=region['total_A_per_um']*c['multiplier']*mult
            err=abs(fd/predicted-1);even=abs((values['plus_'+label]+values['minus_'+label])/2-zero)/abs(fd)
            signal=abs(fd)/max(drift,1e-300)
            signs=(values['plus_'+label]-zero)*predicted>0 and (values['minus_'+label]-zero)*predicted<0
            passed=len(group)==4 and all(r['qualified']=='True' for r in group) and err<=.001 and even<=.01 and linearity<=.001 and signal>=100 and signs and abs(math.log10(zero/float(c['current_A_per_um'])))<=1e-5
            cal.append({'case':c['case'],'amplitude':label,'predicted_A_per_um':predicted,'actual_A_per_um':fd,'relative_error':err,
                'even_fraction':even,'two_amplitude_relative':linearity,'signal_to_zero_drift':signal,
                'response_over_injected_current':fd/(float(c['current_A_per_um'])*c['relative_amplitude']*mult),'passed':passed})
    a.write_csv(OUT/'calibration.csv',cal)
    a.write(OUT/'result.json',{'strict_states':sum(r['qualified']=='True' for r in rows),'passed_amplitudes':sum(r['passed'] for r in cal),'total_amplitudes':4,'m82_released':False,'m83_released':False})
    print(cal,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=('prepare','run','analyze'))
    globals()[parser.parse_args().action]()
