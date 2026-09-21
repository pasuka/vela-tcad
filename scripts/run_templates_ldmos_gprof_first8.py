"""Profile R9 Vg8 first eight exact points with isolated Release controls.

Input bundle must be produced by the hash-checking production exporter. No
physics or convergence setting is changed. Profiles include a cold sweep and
a drain-only replay from the control's initialized zero-drain state.
"""
import argparse
import copy
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

from analyze_templates_ldmos_contact_sweep import numerical
from analyze_templates_ldmos_predictor_study import state_delta
from run_bv_performance_benchmarks import generate_gprof
from run_templates_ldmos_contact_sweep import read, sha
from run_templates_ldmos_electrothermal_curve import state_gate

STATE_KEYS = ('state_interleaved', 'referenced_state_interleaved',
              'electron_qf_reference_V', 'hole_qf_reference_V')


def cost(results):
    total = dict(points=len(results), updates=0, trials=0, performance={})
    for data in results:
        total['updates'] += data['newton_updates']
        total['trials'] += sum(h['line_search_trials'] for h in data['history'])
        perf = data['performance']
        assert perf['linear_solver'] == 'umfpack'
        for key, value in perf.items():
            if isinstance(value, (float, int)):
                total['performance'][key] = total['performance'].get(key, 0) + value
            elif isinstance(value, list) and all(isinstance(v, (int, float)) for v in value):
                old = total['performance'].setdefault(key, [0]*len(value))
                total['performance'][key] = [a+b for a, b in zip(old, value)]
    return total


def inspect(directory, expected):
    ledger = read(directory/'results/ledger.json')
    assert ledger['status'] == 'complete'
    assert [p['bias_V'] for p in ledger['exact_points']] == expected
    for p in ledger['exact_points']:
        assert state_gate(read(p['result']), p['bias_V'])['pass_gate']
    drain = [read(Path(r['directory'])/'output.json') for r in ledger['runs']]
    for r, data in zip(ledger['runs'], drain):
        assert state_gate(data, r['bias_V'])['pass_gate'] == r['gate']['pass_gate']
    init = [read(r['result']) for r in ledger.get('initialization_runs', [])]
    answer = dict(drain=cost(drain), initialization=cost(init),
                  failed_attempts=sum(not r['gate']['pass_gate'] for r in ledger['runs']),
                  drain_attempt_wall_seconds=sum(r['wall_seconds'] for r in ledger['runs']),
                  ledger_wall_seconds=ledger['wall_seconds'],
                  ledger_sha256=sha(directory/'results/ledger.json'))
    return ledger, drain, init, answer


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('bundle', 'release', 'instrumented', 'gprof', 'output'):
        p.add_argument('--'+key, type=Path, required=True)
    a = p.parse_args()
    bundle, out = a.bundle.resolve(), a.output.resolve()
    manifest = read(bundle/'manifest.json')
    for name, digest in manifest['files_sha256'].items():
        assert sha(bundle/name) == digest, name
    cfg, deck = read(bundle/'input_vg8.json'), read(bundle/'vg8.json')
    assert cfg['diagnostic_near_steady_contact_consistency']
    for key in ('reuse_neutral_contact_roots', 'diagnostic_neutral_root_newton',
                'diagnostic_local_qf_limiter'):
        assert not cfg.get(key, False), key
    assert cfg['performance_profiling']
    assert deck['initialization']['mode'] == 'neutral_300K'
    assert deck['initialization']['gate_voltage_V'] == 8
    expected = deck['sweep']['bias_points_V'][:8]
    assert len(expected) == 8 and expected[0] == 0
    deck['sweep']['bias_points_V'] = expected
    out.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    env.pop('GMON_OUT_PREFIX', None)
    report = dict(status='running', exact_bias_V=expected, runs=[],
                  bundle_manifest_sha256=sha(bundle/'manifest.json'),
                  scope='R9 Release cold first8 control, gprof cold and drain-only replay; original gates')

    def save():
        (out/'summary.json').write_text(json.dumps(report, indent=2), encoding='utf-8')

    save()
    baseline = None
    try:
        for name, source in (('release_cold', a.release), ('gprof_cold', a.instrumented),
                             ('gprof_drain', a.instrumented)):
            directory = out/name
            directory.mkdir()
            runtime = directory/'runtime'
            runtime.mkdir()
            source = source.resolve()
            exe = runtime/source.name
            shutil.copy2(source, exe)
            for dll in source.parent.glob('*.dll'):
                shutil.copy2(dll, runtime/dll.name)
            frozen = {f.name: sha(f) for f in runtime.iterdir()}
            run_cfg, run_deck = copy.deepcopy(cfg), copy.deepcopy(deck)
            if name == 'gprof_drain':
                initial = read(baseline[0]['initialized_result'])
                for key in STATE_KEYS:
                    run_cfg[key] = initial[key]
                run_deck['initialization'] = dict(mode='provided_state')
            run_deck['input_file'] = str(directory/'input.json')
            run_deck['output_directory'] = str(directory/'results')
            (directory/'input.json').write_text(json.dumps(run_cfg), encoding='utf-8')
            (directory/'deck.json').write_text(json.dumps(run_deck, indent=2), encoding='utf-8')
            start = time.perf_counter()
            with (directory/'run.log').open('x', encoding='utf-8') as log:
                proc = subprocess.Popen([str(exe), '--config', str(directory/'deck.json')],
                                        cwd=directory, env=env, stdout=log, stderr=subprocess.STDOUT)
                while True:
                    try:
                        code = proc.wait(timeout=30)
                        break
                    except subprocess.TimeoutExpired:
                        progress = dict(run=name, elapsed_seconds=time.perf_counter()-start)
                        path = directory/'results/ledger.json'
                        if path.exists():
                            try:
                                ledger = read(path)
                                progress.update(bias_V=ledger['accepted_bias_V'],
                                                exact_points=len(ledger['exact_points']),
                                                initialization_stages=len(ledger.get('initialization_runs', [])))
                            except (OSError, ValueError):
                                pass
                        print(json.dumps(progress), flush=True)
            wall = time.perf_counter()-start
            assert code == 0, (name, code)
            ledger, drain, init, row = inspect(directory, expected)
            row.update(name=name, external_wall_seconds=wall, runtime_sha256=frozen,
                       input_sha256=sha(directory/'input.json'), deck_sha256=sha(directory/'deck.json'))
            if baseline is None:
                baseline = ledger, drain, init
            else:
                assert len(drain) == len(baseline[1])
                delta = [0.]*4
                for r, old, data, prior in zip(ledger['runs'], baseline[0]['runs'], drain, baseline[1]):
                    assert r['bias_V'] == old['bias_V'] and r['parent_bias_V'] == old['parent_bias_V']
                    delta = [max(x, y) for x, y in zip(delta, state_delta(data, prior))]
                assert all(v <= t for v, t in zip(delta, (1e-8, 1e-8, 1e-8, 1e-7))), delta
                row['drain_max_state_delta_to_release'] = delta
                row['drain_numerical_exact_to_release'] = numerical(drain) == numerical(baseline[1])
                if name == 'gprof_cold':
                    assert len(init) == len(baseline[2])
                    row['initialization_numerical_exact_to_release'] = numerical(init) == numerical(baseline[2])
            if name.startswith('gprof'):
                generate_gprof(directory, exe, a.gprof.resolve())
                row['gmon_sha256'] = sha(directory/'gmon.out')
            for f, digest in frozen.items():
                assert sha(runtime/f) == digest
            report['runs'].append(row)
            save()
            print(json.dumps(dict(run=name, wall_seconds=wall, drain=row['drain']['updates'],
                                  initialization=row['initialization']['updates'])), flush=True)
        report['status'] = 'complete'
        save()
    except BaseException as error:
        report.update(status='failed', error=repr(error))
        save()
        raise


if __name__ == '__main__':
    main()
