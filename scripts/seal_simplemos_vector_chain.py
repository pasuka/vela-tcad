"""Review and seal the isolated vector-drive chain validation."""
import argparse
import math
from pathlib import Path
import validate_simplemos_vector_chain as run

EVIDENCE=run.OUT/'evidence.json'
DOC=run.REPO/'docs/validation/simplemos_vector_chain_validation_2026-09-05.md'


def seal():
    run.p.audit.verify(run.FREEZE)
    fixed=run.p.audit.read(run.OUT/'fixed_state_result.json')
    result=run.p.audit.read(run.OUT/'result.json')
    assert fixed['passed'] and result['passed'] and result['strict_states']==5
    assert not result['m82_released'] and not result['m83_released']
    checks=run.p.audit.rows(run.OUT/'jvp_checks.csv')
    assert len(checks)==12 and all(r['pass']=='True' for r in checks)
    assert len({(r['case'],r['mode'],r['step']) for r in checks})==12
    for mode in ('off','on'):
        for case in ('1e13','5e12'):
            rows=run.p.audit.rows(run.LOCAL/mode/case/'jvp.csv')
            assert len(rows)==53352
            assert len({(r['mode'],r['step'],r['row_block'],r['node_id']) for r in rows})==len(rows)
            assert all(math.isfinite(float(r[k])) for r in rows for k in ('analytic','fd','actual_endpoint_fd'))
    rows=run.p.audit.rows(run.OUT/'dc_ledger.csv')
    assert len(rows)==5 and all(r['qualified']=='True' and int(r['exit_code'])==0 and int(r['local_violations'])==0 for r in rows)
    assert max(abs(float(r['current_relative_change'])) for r in rows)<=1e-8
    for name,_ in run.p.CASES:
        original=run.p.audit.read(run.p.LOCAL/name/'config.json')
        candidate=run.p.audit.read(run.LOCAL/'dc'/name/'config.json')
        original.pop('output_state_file');candidate.pop('output_state_file')
        assert original==candidate,'Changed DC inputs other than output location'
    assert max(float(r['relative_error']) for r in run.p.audit.rows(run.OUT/'dc_calibration.csv'))<=.001
    paths=[Path(__file__),DOC,run.FREEZE]
    paths += [p for d in (run.LOCAL,run.OUT) for p in d.rglob('*') if p.is_file() and p!=EVIDENCE]
    paths += [run.REPO/p for p in run.p.audit.read(run.FREEZE)['input_hashes']]
    run.p.audit.write(EVIDENCE,{'status':'isolated_vector_chain_validated_n23_low_Vd_only','new_nonlinear_solves':5,
        'read_only_Jv_calls':4,'adjoint_calls':1,'read_only_acceptance_calls':5,'new_sentaurus_runs':0,
        'production_changes':False,'m82_released':False,'m83_released':False,
        'input_hashes':{(run.p.audit.rel(p) if p.is_relative_to(run.REPO) else str(p)):run.p.audit.sha(p) for p in sorted(set(paths))}})
    print('Sealed residual/Jv controls, corrected adjoint, and five strict DC states')


def verify():
    run.p.audit.verify(run.FREEZE);run.p.audit.verify(EVIDENCE)
    print('Vector-drive chain validation evidence verified')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=('seal','verify'))
    {'seal':seal,'verify':verify}[parser.parse_args().action]()
