"""Check PTC defect decisions against frozen logs and retain unsuccessful costs."""
import argparse
import json
import math
from pathlib import Path
from analyze_templates_ldmos_pseudo_transient import checked_matrix, sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--matrix', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    summary, runs = checked_matrix(args.matrix)
    expected = {(v, g) for v in ('baseline', 'ptc_v1_s1', 'ptc_defect_ser_s1', 'ptc_defect_model_s1') for g in (4, 8)}
    if len(runs) != 8 or {(r['variant'], r['vg']) for r, _ in runs} != expected:
        raise ValueError('Incomplete acceptance controls')
    result = dict(schema='vela.pseudo_acceptance_findings.v1', matrix_sha256=sha(args.matrix/'summary.json'), runs=[])
    for run, data in runs:
        steps = data.get('pseudo_transient', {}).get('steps', [])
        trials = [t for s in steps for t in s.get('defect_trials', [])]
        for step in steps:
            for t in step.get('defect_trials', []):
                bound = (1.-1e-4*t['alpha'])*step['fixed_residual_before']
                if not all(math.isfinite(t[k]) for k in ('alpha', 'defect_norm', 'model_error', 'steady_fixed_norm')):
                    raise ValueError('Nonfinite defect metric')
                if t['accepted'] != (t['defect_norm'] <= bound):
                    raise ValueError('Armijo decision differs from recorded actual defect')
        ptc = data.get('pseudo_transient', {})
        result['runs'].append(dict(variant=run['variant'], vg=run['vg'], qualified=run['qualified'],
            updates=run['newton_updates'], attempts=len(data['history']), stop=run['stop'],
            trials=run['trials'], assemblies=run['performance']['assembly_calls'],
            factorizations=run['performance']['factorizations'], fallbacks=run['fallbacks'],
            wall_seconds=run['wall_seconds'], defect_evaluations=ptc.get('defect_evaluations', 0),
            defect_seconds=ptc.get('defect_seconds', 0.), defect_trial_checks=len(trials),
            accepted_defect_steps=sum(s.get('defect_used', False) for s in steps),
            accepted_steady_increases=sum(s.get('defect_used', False) and s['fixed_residual_after']>s['fixed_residual_before'] for s in steps),
            initial_tau_s=ptc.get('initial_tau_s'), final_tau_s=ptc.get('next_tau_s'),
            fixed_residual_final_ratio=(steps[-1]['fixed_residual_after']/steps[0]['fixed_residual_before'] if steps else None),
            terminal_blocks=data['electrical_block_gates'], terminal_row_gate=data['carrier_row_gate'],
            max_state_delta=run['max_state_delta']))
    result['scope'] = 'Intermediate defect acceptance is not steady convergence. All failed attempts remain counted.'
    with args.output.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
    print(json.dumps([dict(variant=r['variant'], vg=r['vg'], updates=r['updates'], qualified=r['qualified'],
                          accepted_defect_steps=r['accepted_defect_steps']) for r in result['runs']]))


if __name__ == '__main__':
    main()
