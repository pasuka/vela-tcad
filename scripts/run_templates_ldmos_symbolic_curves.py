"""Serial cold R10/candidate curves; preserve original gates and all failed attempts."""
import argparse
import copy
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

from analyze_templates_ldmos_predictor_study import state_delta
from run_templates_ldmos_contact_sweep import read, sha
from run_templates_ldmos_gprof_first8 import inspect
from run_templates_ldmos_symbolic_study import FLAGS, FULL, REMOTE


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('bundle', 'evidence', 'runner', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--points', type=int, choices=(8, 31), default=8)
    parser.add_argument('--gates', type=int, nargs='+', choices=(4, 8), default=[4, 8])
    parser.add_argument('--variants', nargs='+', choices=tuple(FLAGS), default=['baseline', 'high', 'combined'])
    args = parser.parse_args()
    bundle, root, out = args.bundle.resolve(), args.evidence.resolve(), args.output.resolve()
    manifest = read(bundle/'manifest.json')
    for name, digest in manifest['files_sha256'].items():
        assert sha(bundle/name) == digest, name
    out.mkdir(parents=True, exist_ok=False)
    runtime = out/'runtime'
    runtime.mkdir()
    binary = runtime/args.runner.name
    shutil.copy2(args.runner, binary)
    if os.name == 'nt':
        import ctypes
        ctypes.windll.kernel32.SetErrorMode(0x0001 | 0x0002 | 0x8000)
        for dll in args.runner.parent.glob('*.dll'):
            shutil.copy2(dll, runtime/dll.name)
    frozen = {p.name: sha(p) for p in runtime.iterdir()}
    report = dict(status='running', points=args.points, runtime_sha256=frozen,
                  bundle_manifest_sha256=sha(bundle/'manifest.json'), runs=[])

    def save():
        temporary = out/'summary.tmp'
        temporary.write_text(json.dumps(report, indent=2), encoding='utf-8')
        os.replace(temporary, out/'summary.json')

    def local(path):
        assert path.startswith(REMOTE+'/'), path
        return root/path[len(REMOTE)+1:]

    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    save()
    try:
        for gate in args.gates:
            reference = read(root/FULL/('r0_vg%d_combined/results/ledger.json' % gate))
            assert reference['status'] == 'complete'
            for variant in args.variants:
                directory = out/('vg%d_%s' % (gate, variant))
                directory.mkdir()
                cfg = read(bundle/('input_vg%d.json' % gate))
                deck = copy.deepcopy(read(bundle/('vg%d.json' % gate)))
                assert deck['reuse_static_preparation'] and cfg['reuse_ialmob_thermal_high_field']
                cfg['diagnostic_ialmob_explicit_high_field'], cfg['diagnostic_ialmob_generated_low_field'] = FLAGS[variant]
                deck['sweep']['bias_points_V'] = deck['sweep']['bias_points_V'][:args.points]
                deck.update(input_file=str(directory/'input.json'), output_directory=str(directory/'results'))
                (directory/'input.json').write_text(json.dumps(cfg), encoding='utf-8')
                (directory/'deck.json').write_text(json.dumps(deck, indent=2), encoding='utf-8')
                started = time.perf_counter()
                with (directory/'run.log').open('x', encoding='utf-8') as log:
                    proc = subprocess.Popen([str(binary), '--config', str(directory/'deck.json')],
                                            cwd=directory, env=env, stdout=log, stderr=subprocess.STDOUT)
                    while True:
                        try:
                            code = proc.wait(timeout=30)
                            break
                        except subprocess.TimeoutExpired:
                            progress = dict(gate=gate, variant=variant, elapsed=time.perf_counter()-started)
                            ledger_path = directory/'results/ledger.json'
                            if ledger_path.exists():
                                ledger = read(ledger_path)
                                progress.update(bias_V=ledger['accepted_bias_V'], exact_points=len(ledger['exact_points']))
                            print(json.dumps(progress), flush=True)
                wall = time.perf_counter()-started
                row = dict(gate_V=gate, variant=variant, exit_code=code, external_wall_seconds=wall,
                           input_sha256=sha(directory/'input.json'), deck_sha256=sha(directory/'deck.json'))
                # Record cost before parsing so a malformed or missing output
                # cannot disappear from the report as a cost-free failure.
                row['qualified'] = False
                report['runs'].append(row)
                save()
                if code == 0:
                    ledger, _, _, cost = inspect(directory, deck['sweep']['bias_points_V'])
                    row.update(cost)
                    differences = []
                    for point in ledger['exact_points']:
                        original = next(p for p in reference['exact_points'] if p['bias_V'] == point['bias_V'])
                        differences.append(dict(bias_V=point['bias_V'], maximum=state_delta(read(Path(point['result'])), read(local(original['result'])))))
                    maximum = [max(r['maximum'][k] for r in differences) for k in range(4)]
                    row.update(exact_state_differences=differences, max_state_delta=maximum,
                               qualified=all(x <= limit for x, limit in zip(maximum, (1e-8, 1e-8, 1e-8, 1e-7))))
                else:
                    ledger_path = directory/'results/ledger.json'
                    row.update(qualified=False, ledger=read(ledger_path) if ledger_path.exists() else None)
                save()
                print(json.dumps({k: v for k, v in row.items() if k not in ('ledger', 'exact_state_differences', 'drain', 'initialization')}), flush=True)
                if not row['qualified']:
                    raise RuntimeError('Curve failed qualification: %d %s' % (gate, variant))
        assert all(sha(runtime/name) == digest for name, digest in frozen.items())
        report.update(status='complete', all_qualified=True)
        save()
    except BaseException as error:
        report.update(status='failed', error=repr(error))
        save()
        raise


if __name__ == '__main__':
    main()
