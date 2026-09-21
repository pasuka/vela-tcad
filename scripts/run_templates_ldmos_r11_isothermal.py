"""Serial G3 and full zero-drain D5/D4 regressions of the frozen R11 source.

Use the existing isothermal physics and gates. The R11 electrothermal-only
kernel switches are not inserted into the separate dc_sweep input schema.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

from run_templates_ldmos_linked_d5 import cpu_times, preflight, read, write
from analyze_templates_ldmos_stage4_d5 import analyze

ROOT = Path(__file__).resolve().parents[1]
G3_PRIOR = ROOT/'reference_staging/templates_ldmos_joint_20260914/regression_r6/g3_auger'
PROFILES = ROOT/'reference_tcad/templates_ldmos_sentaurus2022/profiles'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def relocate(value, old, new):
    if isinstance(value, dict):
        return {k: relocate(v, old, new) for k, v in value.items()}
    if isinstance(value, list):
        return [relocate(v, old, new) for v in value]
    if isinstance(value, str):
        for source in (str(old), old.as_posix()):
            if value.startswith(source + '/') or value.startswith(source + '\\'):
                return str(new) + value[len(source):]
    return value


def csv_rows(path):
    with Path(path).open(encoding='utf-8', newline='') as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runner', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if os.name != 'nt' or not __debug__:
        raise RuntimeError('This protocol uses Windows process CPU accounting and active assertions')
    out = args.output.resolve()
    assert not out.exists(), 'Refusing to overwrite evidence'
    runner = args.runner.resolve()
    assert runner.is_file()
    flags = (ROOT/'build-release/build.ninja').read_text(encoding='utf-8')
    assert 'FLAGS = -O3 -DNDEBUG -std=c++20' in flags
    assert '-pg ' not in flags
    # The latest G3 input and its prior frozen dependencies are checked before
    # changing output paths. Neither old exact-point results nor scripts run here.
    prior_plan = read(G3_PRIOR/'plan.json')
    for name, digest in prior_plan['frozen_files'].items():
        assert sha(name) == digest, name
    prior_cfg = read(G3_PRIOR/'control.json')
    assert prior_cfg['solver']['auger_with_generation'] is False
    assert len(prior_cfg['sweep']['bias_points']) == 31
    g3_reference = Path(read(G3_PRIOR/'qualification.json')['reference'])
    bundles = {p: PROFILES/('linked_'+p.lower()+'_auger_no_generation_inputs.json') for p in ('D5', 'D4')}
    for bundle in bundles.values():
        for name, digest in read(bundle)['files'].items():
            assert sha(ROOT/name) == digest, name
    assert shutil.disk_usage(out.parent).free > 12*1024**3, 'Insufficient evidence disk space'
    out.mkdir()
    binary = out/'binary'
    binary.mkdir()
    frozen_runner = binary/runner.name
    shutil.copy2(runner, frozen_runner)
    sources = subprocess.check_output(['git', '-c', 'core.fsmonitor=false', 'ls-files',
        'src', 'include', 'CMakeLists.txt'], cwd=ROOT, text=True).splitlines()
    sources += ['scripts/'+p.name for p in (ROOT/'scripts').glob('*.py')]
    frozen = {}
    for name in sorted(set(sources)):
        dst = out/'source'/name
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT/name, dst)
        frozen[str(dst)] = sha(dst)
    for name in ('CMakeCache.txt', 'build.ninja'):
        shutil.copy2(ROOT/'build-release'/name, binary/name)
    manifest = dict(runner=str(frozen_runner), runner_sha256=sha(frozen_runner),
        backend='Eigen SparseLU/COLAMD (VELA_LINEAR_SOLVER=sparselu)',
        source_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        build_type='Release -O3 -DNDEBUG', frozen_sources=frozen,
        runtime_sha256={str(p): sha(p) for p in (Path('D:/msys64/ucrt64/bin')/n for n in
            ('libwinpthread-1.dll', 'libgcc_s_seh-1.dll', 'libstdc++-6.dll'))})
    write(binary/'manifest.json', manifest)
    for bundle in bundles.values():
        for gate in (4, 8):
            preflight(bundle, ROOT, binary/'manifest.json', gate)
    env = dict(os.environ, VELA_LINEAR_SOLVER='sparselu', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    env['PATH'] = 'D:/msys64/ucrt64/bin;D:/msys64/usr/bin;' + env.get('PATH', '')
    env.pop('GMON_OUT_PREFIX', None)
    report = dict(status='running', source_commit=manifest['source_commit'],
        runner_sha256=manifest['runner_sha256'], backend=manifest['backend'], cases=[],
        scope='G3 prebiased full 31-point sweep; D5/D4 zero-drain seed requalification then full 31-point continuation; original gates',
        inputs_sha256={str(p): sha(p) for p in [G3_PRIOR/'control.json', g3_reference, *bundles.values()]})

    def save():
        write(out/'summary.json', report)

    def execute(name, argv, directory):
        assert sha(frozen_runner) == manifest['runner_sha256']
        start = time.perf_counter()
        row = dict(name=name, argv=list(map(str, argv)), status='running')
        report['cases'].append(row)
        report['active_case'] = name
        with (out/(name+'.log')).open('x', encoding='utf-8') as log:
            proc = subprocess.Popen(list(map(str, argv)), cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
            row['pid'] = proc.pid
            save()
            while proc.poll() is None:
                try:
                    proc.wait(timeout=30)
                except subprocess.TimeoutExpired:
                    row['elapsed_seconds'] = time.perf_counter()-start
                    save()
                    print(json.dumps(dict(active_case=name, elapsed_seconds=row['elapsed_seconds'])), flush=True)
            row.update(exit_code=proc.returncode, process_cpu_seconds=cpu_times(proc),
                wall_seconds=time.perf_counter()-start, status='completed' if proc.returncode == 0 else 'failed')
        save()
        assert proc.returncode == 0, (name, proc.returncode)
        return row

    try:
        save()
        g3 = out/'g3'
        g3.mkdir()
        cfg = relocate(prior_cfg, G3_PRIOR, g3)
        assert relocate(cfg, g3, G3_PRIOR) == prior_cfg
        write(g3/'control.json', cfg)
        row = execute('g3', [frozen_runner, '--config', g3/'control.json'], g3)
        assert len(csv_rows(cfg['output_csv'])) == 31
        assert all(x['converged'] == '1' for x in csv_rows(cfg['output_csv']))
        with (g3/'score.log').open('x', encoding='utf-8') as log:
            subprocess.check_call([sys.executable, str(ROOT/'scripts/analyze_templates_ldmos_g3_idvg.py'),
                '--reference', str(g3_reference), '--candidate', cfg['output_csv'],
                '--balance', cfg['sweep']['diagnostics']['srh_balance']['csv_file'],
                '--output-json', str(g3/'qualification.json'), '--output-md', str(g3/'qualification.md')],
                cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
        qualification = read(g3/'qualification.json')
        assert qualification['status'] == 'pass' and qualification['metrics']['exact_shared_points'] == 31
        row['qualification'] = qualification
        save()
        for profile in ('D5', 'D4'):
            bundle = read(bundles[profile])
            for gate in (4, 8):
                name = profile.lower()+'_vg'+str(gate)
                dest = out/name
                row = execute(name, [sys.executable, ROOT/'scripts/run_templates_ldmos_linked_d5.py',
                    '--physics-profile', profile, '--bundle', bundles[profile], '--manifest', binary/'manifest.json',
                    '--output', dest, '--gate', gate, '--points', '31'], dest)
                assert read(dest/'progress.json')['status'] == 'completed'
                row['audit'] = read(dest/'audit_summary.json')
                assert row['audit']['integrity_pass'] and row['audit']['exact_points'] == 31
                save()
            qualification = analyze(
                {'Vg'+str(g): ROOT/bundle['references'][str(g)] for g in (4, 8)},
                {'Vg'+str(g): out/(profile.lower()+'_vg'+str(g))/'score/curve.csv' for g in (4, 8)},
                {'Vg'+str(g): out/(profile.lower()+'_vg'+str(g))/'score/terminal_balance.csv' for g in (4, 8)},
                out/(profile.lower()+'_joint'), physics_profile=profile)
            assert qualification['engineering']['status'] == qualification['final']['status'] == 'pass'
            report[profile+'_qualification'] = qualification
            save()
        for name, digest in frozen.items():
            assert sha(name) == digest
        report.update(status='pass', exact_points=155)
        report.pop('active_case', None)
    except BaseException as error:
        report.update(status='failed', error=str(error))
        raise
    finally:
        save()
    print('R11_ISOTHERMAL_COMPLETE', flush=True)


if __name__ == '__main__':
    main()
