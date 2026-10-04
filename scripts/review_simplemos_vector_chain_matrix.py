"""Review raw four-point evidence and preserve the physical Id discrepancy ledger."""
import math
from pathlib import Path
import validate_simplemos_vector_chain_matrix as m

a = m.a


def main():
    a.verify(m.FREEZE)
    m.sealed.verify()
    cfg = a.read(m.CONTRACT)
    cal = a.rows(m.OUT / 'calibration.csv')
    assert len(cal) == 16
    assert len({(r['case'], r['mode'], r['amplitude_cm_3']) for r in cal}) == 16
    jvp = a.rows(m.OUT / 'jvp_checks.csv')
    assert len(jvp) == 18
    assert len({(r['case'], r['amplitude_cm_3'], r['step']) for r in jvp}) == 18
    old_errors = {r['case']: float(r['strict_precision_error_dex']) for r in a.rows(m.spatial.precision.OUT / 'current_error_ledger.csv')}
    physics, iterations, observer = [], [], []
    for c in cfg['cases']:
        key = c['case']
        currents = {}
        for mode in ('off', 'on'):
            root = (m.p.LOCAL if mode == 'off' else m.chain.LOCAL / 'dc') if key == m.REUSED else m.LOCAL / key / mode
            for name, _ in m.p.CASES:
                dest = root / name
                filename = 'solve.status.json' if key == m.REUSED and mode == 'off' else 'config.status.json'
                s = a.read(dest / filename)
                audit = a.read(dest / 'acceptance.status.json')
                assert s['exit_code'] == 0 and s['converged']
                assert audit['exit_code'] == 0
                assert audit['carrier_row_convergence']['satisfied']
                assert audit['global_continuity_closure']['satisfied']
                cc = s['contact_currents_A_per_um']
                assert abs(math.fsum(cc.values())) / abs(cc['drain']) <= 1e-8
                if name == 'zero':
                    currents[mode] = cc['drain']
                iterations.append({'case': key, 'mode': mode, 'name': name, 'iterations': s['iterations'], 'reason': s['convergence_reason']})
                deck = a.read(dest / 'config.json')
                deck.pop('output_state_file')
                if mode == 'off':
                    # Same physical inputs, initial state, solver gates for on/off.
                    other = (m.chain.LOCAL / 'dc' / name) if key == m.REUSED else m.LOCAL / key / 'on' / name
                    comparison = a.read(other / 'config.json')
                    comparison.pop('output_state_file')
                    assert deck == comparison
        if key != m.REUSED:
            status = a.read(m.LOCAL / key / 'adjoint/config.status.json')
            identity = abs(status['current_A_per_um'] / float(c['current_A_per_um']) - 1)
            observer.append({'case': key, 'current_identity_relative': identity, 'adjoint_relative_residual': status['adjoint_relative_residual'], 'passed': identity <= 1e-8 and status['adjoint_relative_residual'] <= 1e-10})
            for suffix in ('1e13', '5e12'):
                dest = m.LOCAL / key / 'jvp' / suffix
                a.verify(dest / 'freeze.json')
                for mode in ('off', 'on'):
                    rows = a.rows(dest / mode / 'jvp.csv')
                    assert len(rows) == c['nodes'] * 36
                    assert len({(r['mode'], r['step'], r['row_block'], r['node_id']) for r in rows}) == len(rows)
                    assert all(math.isfinite(float(r[k])) for r in rows for k in ('analytic', 'fd', 'actual_endpoint_fd'))
        for mode in ('off', 'on'):
            err = old_errors[key] + math.log10(currents[mode] / float(c['current_A_per_um']))
            response = next(r for r in cal if r['case'] == key and r['mode'] == mode and float(r['amplitude_cm_3']) == 1e13)
            physics.append({'case': key, 'device': c['device'], 'vd': float(c['vd']), 'vg': float(c['vg']), 'mode': mode,
                'current_A_per_um': currents[mode], 'error_dex': err, 'error_percent': 100 * (10**err - 1),
                'normalized_charge_response_at_1e13': float(response['central_response_A_per_um']) / currents[mode]})
    pairs = []
    for vd in (.05, 1.):
        for mode in ('off', 'on'):
            lo = next(r for r in physics if r['device'] == 'n19' and r['vd'] == vd and r['mode'] == mode)
            hi = next(r for r in physics if r['device'] == 'n23' and r['vd'] == vd and r['mode'] == mode)
            pairs.append({'vd': vd, 'mode': mode, 'high_minus_low_error_dex': hi['error_dex'] - lo['error_dex'],
                'high_to_low_normalized_charge_response': hi['normalized_charge_response_at_1e13'] / lo['normalized_charge_response_at_1e13']})
    a.write_csv(m.OUT / 'physical_error_ledger.csv', physics)
    a.write_csv(m.OUT / 'paired_error_ledger.csv', pairs)
    a.write_csv(m.OUT / 'iteration_ledger.csv', iterations)
    a.write_csv(m.OUT / 'observer_checks.csv', observer)
    result = {'reviewed_states': len(iterations), 'new_states': 30, 'reused_states': 10,
        'dc_calibration_passed': a.read(m.OUT / 'result.json')['passed'],
        'jvp_passed': all(r['passed'] == 'True' for r in jvp),
        'observer_passed': all(r['passed'] for r in observer),
        'no_production_changes': True, 'm82_released': False, 'm83_released': False,
        'input_hashes': {a.rel(f): a.sha(f) for f in (Path(__file__).resolve(),
            m.REPO / 'tests/regression/test_simplemos_vector_chain_matrix.py',
            m.spatial.precision.OUT / 'current_error_ledger.csv')}}
    a.write(m.OUT / 'review.json', result)
    print(result, flush=True)


if __name__ == '__main__':
    main()
