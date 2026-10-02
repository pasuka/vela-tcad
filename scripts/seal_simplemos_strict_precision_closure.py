"""Review and seal the strict precision closure and spatial response evidence."""
import argparse
import copy
import math
from pathlib import Path
import run_simplemos_strict_flux_precision as precision
import trace_simplemos_strict_rejection as trace
import compare_simplemos_qualified_charge_response as response

audit=precision.audit
REPO,ROOT=precision.REPO,precision.ROOT
EVIDENCE=ROOT/'simplemos_strict_precision_closure_evidence.json'
DOC=REPO/'docs/validation/simplemos_strict_rejection_and_precision_closure_2026-09-05.md'


def review():
    if EVIDENCE.exists():raise FileExistsError(EVIDENCE)
    trace.verify();audit.verify(precision.FREEZE);response.verify()
    cases=audit.rows(precision.OUT/'case_ledger.csv')
    assert len(cases)==4 and all(c['qualified']=='True' for c in cases)
    for c in cases:
        before=audit.read(precision.previous.LOCAL/c['case']/'config.json')
        after=audit.read(precision.LOCAL/c['case']/'config.json')
        before.pop('output_state_file');after.pop('output_state_file')
        assert after['solver']['bandgap_narrowing'].pop('equal_ni_flux_evaluation')=='compensated_log_expm1'
        assert before==after,'More than the declared precision axis changed'
    baseline=audit.rows(response.fixed.upstream.OUT/'m80_case_ledger.csv')
    by={(r['device'],float(r['drain_voltage_V'])):r for r in baseline}
    error_rows=[];pair_rows=[]
    for c in cases:
        b=by[c['device'],float(c['vd'])];ref=float(b['sentaurus_current_A_per_um'])
        old_error=math.log10(float(b['vela_current_A_per_um'])/ref)
        error_rows.append({'case':c['case'],'device':c['device'],'vd':c['vd'],
            'previous_error_dex':old_error,'strict_precision_error_dex':math.log10(float(c['current_A_per_um'])/ref),
            'observed_change_dex':float(c['delta_log10_Id_dex'])})
    for vd in (.05,1.):
        lo=next(r for r in error_rows if r['device']=='n19' and float(r['vd'])==vd)
        hi=next(r for r in error_rows if r['device']=='n23' and float(r['vd'])==vd)
        pair_rows.append({'vd':vd,'previous_pair_error_dex':hi['previous_error_dex']-lo['previous_error_dex'],
            'qualified_pair_error_dex':hi['strict_precision_error_dex']-lo['strict_precision_error_dex'],
            'qualified_pair':True})
    audit.write_csv(precision.OUT/'current_error_ledger.csv',error_rows)
    audit.write_csv(precision.OUT/'paired_error_ledger.csv',pair_rows)
    normalized=[];calibration=[]
    for r in audit.rows(response.OUT/'response_ledger.csv'):
        if not r['sentaurus_green_delta_Id_A_per_um']:continue
        c=next(c for c in cases if c['case']==r['case']);ref=float(by[r['device'],float(r['vd'])]['sentaurus_current_A_per_um'])
        normalized.append({'support':r['support'],'device':r['device'],'vd':r['vd'],
            'native_dlog10_Id':float(r['sentaurus_green_delta_Id_A_per_um'])/(ref*math.log(10)),
            'vela_dlog10_Id':float(r['vela_delta_log10_Id_linear']),
            'normalized_response_ratio':float(r['vela_to_native_response_ratio'])*ref/float(c['current_A_per_um']),
            'native_direction_fd_qualified':r['native_direction_fd_qualified']})
    for c in cases:
        s=audit.read(response.LOCAL/c['case']/'config.status.json');v=s['electron_volume_response']
        # Reuse M79's frozen signal definition; never classify null hole-current directions as passed FD tests.
        floor=max(1e-20,1e-12*s['state_derivative_norm'])
        good=[r for r in v['finite_difference_checks'] if abs(r['current_analytic'])>floor]
        maxj=max(r['jvp_relative_error'] for r in v['finite_difference_checks'])
        maxg=max(r['current_relative_error'] for r in good)
        assert good and maxj<=1e-3 and maxg<=1e-3 and v['duality_relative_error']<=1e-8 and v['linear_solve_relative_residual']<=1e-8
        calibration.append({'case':c['case'],'jvp_relative_max':maxj,'current_fd_signal_floor':floor,
            'qualified_current_fd_checks':len(good),'unqualified_current_fd_checks':len(v['finite_difference_checks'])-len(good),
            'current_fd_relative_max_qualified':maxg,'duality_relative':v['duality_relative_error'],'linear_solve_relative':v['linear_solve_relative_residual']})
    audit.write_csv(response.OUT/'normalized_response.csv',normalized)
    audit.write_csv(response.OUT/'numerical_calibration.csv',calibration)
    paths=[Path(__file__),trace.SCRIPT,precision.SCRIPT,response.SCRIPT,DOC]
    for folder in (trace.LOCAL,trace.OUT,precision.LOCAL,precision.OUT,response.LOCAL,response.OUT):
        paths += [p for p in folder.rglob('*') if p.is_file()]
    audit.write(EVIDENCE,{'status':'four_strict_states_qualified_spatial_directions_screened',
        'new_nonlinear_reclosures':7,'new_sentaurus_runs':0,'new_full_DD_adjoint_probes':4,
        'strict_qualified':4,'full_curve_qualified':False,'m82_released':False,
        'input_hashes':{audit.rel(p):audit.sha(p) for p in sorted(set(paths))}})
    print('Sealed four qualified states, unchanged acceptance profiles, and limited spatial screening')


def verify():
    trace.verify();audit.verify(precision.FREEZE);response.verify();audit.verify(EVIDENCE)
    print('Strict precision closure evidence verified')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('review','verify'))
    {'review':review,'verify':verify}[p.parse_args().action]()
