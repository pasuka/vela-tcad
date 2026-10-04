"""Read-only cross-check of raw local-flux outputs and frozen acceptance gates."""
import math
from pathlib import Path
import validate_simplemos_local_conservative_flux as v
import resolve_simplemos_contact_flux_signal as extra

a = v.a


def main():
    for manifest in (v.FREEZE, extra.FREEZE, v.OUT/'field_audit_evidence.json'):
        a.verify(manifest)
    states = []
    extractors = []
    responses = []
    for c in a.read(v.CONTRACT)['cases']:
        for region in ('zero', 'drain_contact', 'channel'):
            s = a.read(v.LOCAL/c['case']/'preflight'/region/'functional.status.json')
            error = abs(s['current_A_per_um']-s['contact_current_extractor_A_per_um'])/abs(s['current_A_per_um'])
            assert error < 1e-8
            extractors.append(dict(case=c['case'], region=region, relative_error=error))
        for root, supplemental in ((v.LOCAL, False), (extra.LOCAL, True)):
            configs = sorted((root/c['case']).glob('*/config.json' if supplemental else '**/config.json'))
            assert len(configs) == (4 if supplemental else 9)
            for config in configs:
                dest = config.parent
                s = a.read(dest/'config.status.json')
                audit = a.read(dest/'acceptance.status.json')
                cc = s['contact_currents_A_per_um']
                kcl = abs(math.fsum(cc.values()))/abs(cc['drain'])
                ok = s['exit_code'] == audit['exit_code'] == 0 and s['converged'] and audit['carrier_row_convergence']['satisfied'] and audit['global_continuity_closure']['satisfied'] and kcl <= 1e-8
                assert ok
                states.append(dict(path=a.rel(dest), original_acceptance=ok, kcl_over_Id=kcl,
                    electron_global_qualified=audit['global_continuity_closure']['electron']['qualified'],
                    hole_global_qualified=audit['global_continuity_closure']['hole']['qualified']))
    assert len(states) == 26
    for folder, supplemental in ((v.OUT, False), (extra.OUT, True)):
        contract = a.read(folder/'contract.json')
        for row in a.rows(folder/'calibration.csv'):
            c = next(c for c in contract['cases'] if c['case'] == row['case'])
            name = 'drain_contact' if supplemental else row['region']
            region = next(r for r in c['regions'] if r['name'] == name)
            root = (extra.LOCAL if supplemental else v.LOCAL)/c['case']
            if not supplemental:
                root = root/name
            amp = row['amplitude']
            plus = a.read(root/('plus_'+amp)/'config.status.json')['contact_currents_A_per_um']['drain']
            minus = a.read(root/('minus_'+amp)/'config.status.json')['contact_currents_A_per_um']['drain']
            fd = (plus-minus)/2
            assert fd == float(row['actual_A_per_um'])
            mult = (1 if amp == 'full' else .5)*(c['multiplier'] if supplemental else 1)
            direct = region['direct_A_per_um']*mult
            feedback = region['feedback_A_per_um']*mult
            zero = a.read(v.LOCAL/c['case']/'zero/config.status.json')['contact_currents_A_per_um']['drain']
            failures = []
            for key, limit in (('relative_error', .001), ('two_amplitude_relative', .001), ('even_fraction', .01)):
                if float(row[key]) > limit:
                    failures.append(key)
            if float(row['signal_to_zero_drift']) < 100:
                failures.append('signal_to_zero_drift')
            if not ((plus-zero)*float(row['predicted_A_per_um']) > 0 and (minus-zero)*float(row['predicted_A_per_um']) < 0):
                failures.append('signed_response')
            if supplemental and any(r['qualified'] != 'True' for r in a.rows(folder/'dc.csv') if r['case'] == c['case']):
                failures.append('additional_state_guard')
            assert (not failures) == (row['passed'] == 'True')
            responses.append(dict(case=c['case'], supplemental=supplemental, region=name, amplitude=amp,
                passed=row['passed'] == 'True', failures=failures,
                direct_A_per_um=direct, feedback_prediction_A_per_um=feedback,
                feedback_fd_A_per_um=fd-direct,
                feedback_only_relative_error=abs((fd-direct)/feedback-1),
                feedback_only_does_not_qualify_net_response=True))
    assert sum(r['passed'] for r in responses) == 4
    assert all(r['region'] == 'channel' for r in responses if r['passed'])
    a.write(v.OUT/'review.json', dict(input_hashes={a.rel(Path(__file__).resolve()):a.sha(Path(__file__).resolve())},
        original_acceptance_passed=26, states=states, extractor_checks=extractors, responses=responses,
        supplemental_state_guard_passed=6, response_amplitudes_passed=4, response_amplitudes_total=12,
        qualification='Only the prescribed Vela channel edge continuity perturbation is qualified at both biases.',
        new_sentaurus_runs=0, m82_released=False, m83_released=False))
    print('Raw-output review passed: 26 original acceptance checks, 6 extractor checks; 4/12 response amplitudes qualified.')


if __name__ == '__main__':
    main()
