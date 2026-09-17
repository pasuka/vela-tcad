"""Audit isolated neutral-root Newton representative pairs and repeatability."""
import argparse
import hashlib
import json
from pathlib import Path

from analyze_templates_ldmos_contact_sweep import numerical
from run_templates_ldmos_contact_sweep import read, sha
from run_templates_ldmos_electrothermal_curve import state_gate


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--matrix', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--repeats', type=int, default=2)
    a = p.parse_args()
    if a.repeats < 1:
        p.error('Positive repeat count required')
    r = read(a.matrix/'summary.json')
    variants = ('baseline', 'neutral_newton')
    expected = {(i, c, v) for i in range(a.repeats) for c in range(8) for v in variants}
    assert r['status'] == 'complete' and len(r['runs']) == len(expected)
    assert {(x['repeat'], x['case'], x['variant']) for x in r['runs']} == expected
    for name, digest in r['runtime_sha256'].items():
        assert sha(a.matrix/name) == digest, name
    first = {}; signatures = {}; inputs = {}; hashes = {}; delta = [0.]*4
    for x in r['runs']:
        key = x['repeat'], x['case'], x['variant']
        directory = a.matrix/('r%d_case%d_%s' % key)
        assert sha(directory/'input.json') == x['input_sha256']
        digest = sha(directory/'output.json')
        assert digest == x['output_sha256']; hashes[str(directory/'output.json')] = digest
        data = read(directory/'output.json'); inputs[key] = read(directory/'input.json')
        assert x['qualified'] and x['numerical_exact_to_baseline']
        assert state_gate(data, x['bias_V'])['pass_gate']
        delta = [max(v, w) for v, w in zip(delta, x['max_state_delta'])]
        signature = hashlib.sha256(json.dumps(numerical(data), sort_keys=True,
                                             separators=(',', ':')).encode()).hexdigest()
        signatures[key] = signature
        case = x['case'], x['variant']
        if x['repeat'] == 0:
            first[case] = signature
        else:
            assert first[case] == signature, 'Nonrepeatable numerical trajectory'
    assert all(v <= t for v, t in zip(delta, (1e-8, 1e-8, 1e-8, 1e-7)))
    for i in range(a.repeats):
        for c in range(8):
            original = inputs[i, c, 'baseline']
            assert not original.get('reuse_neutral_contact_roots', False)
            assert not original.get('diagnostic_neutral_root_newton', False)
            assert inputs[i, c, 'neutral_newton'] == dict(original, diagnostic_neutral_root_newton=True)
            assert signatures[i, c, 'baseline'] == signatures[i, c, 'neutral_newton']
    report = dict(status='pass', all_original_gates=True, exact_to_baseline=True,
                  repeats=a.repeats, repeated_numerics_exact=a.repeats > 1,
                  max_state_delta_to_frozen_reference=delta, groups=[], hashes=hashes)
    for i in range(a.repeats):
        for v in variants:
            xs = [x for x in r['runs'] if x['repeat'] == i and x['variant'] == v]
            q = dict(repeat=i, variant=v, updates=sum(x['updates'] for x in xs),
                     trials=sum(x['trials'] for x in xs),
                     assemblies=sum(x['performance']['assembly_calls'] for x in xs),
                     wall_seconds=sum(x['wall_seconds'] for x in xs),
                     assembly_seconds=sum(x['performance']['assembly_seconds'] for x in xs),
                     factorization_seconds=sum(x['performance']['factorization_seconds'] for x in xs),
                     root_solves=sum(x['performance']['neutral_root_counts'][0] for x in xs),
                     root_iteration_counts=[sum(x['performance']['neutral_root_iteration_counts'][k]
                                                for x in xs) for k in range(5)])
            report['groups'].append(q)
    a.output.write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps({k:v for k, v in report.items() if k != 'hashes'}, indent=2))


if __name__ == '__main__':
    main()
