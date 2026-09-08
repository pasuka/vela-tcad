"""Seal successful DC states and the explicitly failed adjoint FD qualification."""
import argparse
import math
from pathlib import Path
import validate_simplemos_vela_channel_charge as run

EVIDENCE = run.OUT / 'evidence.json'
DOC = run.REPO / 'docs/validation/simplemos_vela_channel_charge_fd_execution_2026-09-05.md'


def review():
    run.audit.verify(run.FREEZE)
    cfg = run.audit.read(run.CONTRACT)
    cases = run.audit.rows(run.OUT/'case_ledger.csv')
    cal = run.audit.rows(run.OUT/'calibration.csv')
    result = run.audit.read(run.OUT/'result.json')
    assert len(cases)==5 and len(cal)==2
    assert {r['name'] for r in cases}=={n for n,_ in run.CASES}
    assert result['qualified_states']==sum(r['qualified']=='True' for r in cases)==5
    assert result['passed_amplitudes']==sum(r['pass']=='True' for r in cal)==0
    assert all(float(r['relative_error'])>cfg['gates']['fd_vs_adjoint_relative'] for r in cal)
    assert result['status']=='failed' and not result['m82_released'] and not result['m83_released']
    for name,_ in run.CASES:
        deck = run.audit.read(run.LOCAL/name/'config.json')
        assert deck['solver']==cfg['solver']
        state = run.audit.rows(run.LOCAL/name/'state.csv')
        assert len(state)==len({r['node_id'] for r in state})
        assert all(math.isfinite(float(v)) for r in state for v in r.values())
    source = run.audit.rows(run.OUT/'source_preflight.csv')
    assert len(source)==5 and all(float(r['source_relative_L2'])<=cfg['gates']['source_relative_L2'] for r in source)
    assert all(int(r['nonzero_nodes'])==(0 if r['name']=='zero' else 186) for r in source)
    field = run.audit.rows(run.OUT/'field_response_nodes.csv')
    assert len(field)==1884
    assert len({(r['amplitude_cm_3'],r['node_id']) for r in field})==1884
    assert len(run.audit.rows(run.OUT/'field_response_summary.csv'))==18
    assert len(run.audit.rows(run.OUT/'actual_direction_gradient.csv'))==2
    assert len(run.audit.rows(run.OUT/'stationarity_projection.csv'))==2
    print('Reviewed: 5 strict states; 0/2 adjoint response calibrations; source and field coverage intact')


def seal():
    review()
    paths = [Path(__file__), DOC, run.FREEZE]
    for folder in (run.LOCAL,run.OUT):
        paths += [p for p in folder.rglob('*') if p.is_file() and p!=EVIDENCE]
    paths += [run.REPO/'scripts'/n for n in ('compare_simplemos_channel_fd_fields.py',
        'audit_simplemos_channel_response_gradient.py','audit_simplemos_channel_stationarity.py')]
    paths += [run.REPO/p for p in run.audit.read(run.FREEZE)['input_hashes']]
    run.audit.write(EVIDENCE,{'status':'five_strict_states_passed_adjoint_direction_calibration_failed',
        'new_nonlinear_solves':5,'new_sentaurus_runs':0,'new_bias_sweeps':0,
        'source_residual_probes':6,'read_only_acceptance_probes':5,'final_residual_probes':5,
        'm82_released':False,'m83_released':False,
        'input_hashes':{(run.audit.rel(p) if p.is_relative_to(run.REPO) else str(p)):run.audit.sha(p) for p in sorted(set(paths))}})
    print('Sealed DC, source, gradient, stationarity and field evidence')


def verify():
    run.audit.verify(run.FREEZE)
    run.audit.verify(EVIDENCE)
    print('Vela channel charge FD evidence verified')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('review','seal','verify'))
    {'review':review,'seal':seal,'verify':verify}[parser.parse_args().action]()
