"""Seal actual-direction defect measurements and the unchanged switch control."""
import argparse
import math
from pathlib import Path
import check_simplemos_actual_direction_jvp as run

EVIDENCE=run.OUT/'evidence.json'


def seal():
    run.prior.audit.verify(run.FREEZE)
    switch=run.OUT/'field_switch_contract.json'
    run.prior.audit.verify(switch)
    result=run.prior.audit.read(run.OUT/'result.json')
    assert result['defect_measurement_qualified'] and len(result['checks'])==2
    cfg=run.prior.audit.read(run.CONTRACT)
    for r in result['checks']:
        assert r['pass'] and all(math.isfinite(r[k]) and r[k]<=v for k,v in cfg['gates'].items())
    control=run.prior.audit.read(run.OUT/'field_switch_result.json')
    assert control['rows']==53352 and control['residual_differences_identical'] and control['Jv_identical']
    rows=run.prior.audit.rows(run.OUT/'row_ledger.csv')
    assert len(rows)==106704
    assert len({(r['case'],r['mode'],r['step'],r['row_block'],r['node_id']) for r in rows})==len(rows)
    for r in rows:
        assert all(math.isfinite(float(r[k])) for k in ('analytic','fd','actual_endpoint_fd','lambda','weighted_defect_A_per_um'))
    assert not result['m82_released'] and not result['m83_released']
    paths=[Path(__file__),Path(run.__file__),run.DOC,run.FREEZE,switch]
    paths += [p for d in (run.LOCAL,run.OUT) for p in d.rglob('*') if p.is_file() and p!=EVIDENCE]
    for contract in (run.FREEZE,switch):
        paths += [run.REPO/p for p in run.prior.audit.read(contract)['input_hashes']]
    run.prior.audit.write(EVIDENCE,{'status':'actual_direction_Jv_defect_localized_electron_phin_block',
        'defect_measurement_qualified':True,'read_only_runner_calls':3,'new_nonlinear_solves':0,'new_sentaurus_runs':0,
        'm82_released':False,'m83_released':False,
        'input_hashes':{(run.prior.audit.rel(p) if p.is_relative_to(run.REPO) else str(p)):run.prior.audit.sha(p) for p in sorted(set(paths))}})
    print('Sealed 106704 primary Jv rows and 53352 switch control rows')


def verify():
    for p in (run.FREEZE,run.OUT/'field_switch_contract.json',EVIDENCE): run.prior.audit.verify(p)
    print('Actual-direction Jv and switch control evidence verified')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('seal','verify'))
    {'seal':seal,'verify':verify}[p.parse_args().action]()
