"""Serial D0 Release closeout: short, full, native fields, restarts, repeats.

Uses a new output directory. Input references are read-only; a failure stops
subsequent stages. The default case omits all linear-policy/root-method keys.
The control disables only cross-point linear-object reuse, not within-point
symbolic reuse. Timings include initialization, rejected attempts and output.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import traceback

from analyze_templates_ldmos_d0 import analyze
from electrothermal_state import read_bound_record
from run_templates_ldmos_electrothermal_curve import state_gate
from run_templates_ldmos_screening_candidates import read, save, digest
from run_templates_ldmos_screening_repeats import compare_states

COUNTERS = ('assembly_seconds', 'assembly_calls', 'factorization_seconds',
            'factorizations', 'symbolic_analyses', 'linear_solve_seconds',
            'ialmob_screening_candidate_calls', 'ialmob_screening_function_evaluations',
            'ialmob_screening_fallbacks', 'linear_object_reused')
POLICY_KEYS = ('electrothermal_linear_solver', 'reuse_linear_analysis',
               'reuse_sparselu_symbolic', 'diagnostic_ialmob_screening_method')


def policy_input(original, enabled):
    result = json.loads(json.dumps(original))
    for key in POLICY_KEYS:
        result.pop(key, None)
    result.get('mobility_SI', {}).get('ialmob', {}).pop('screening_method', None)
    if not enabled:
        result['reuse_linear_analysis'] = False
    return result


def inspect(folder, points, enabled):
    ledger = read(folder/'results/ledger.json')
    if ledger['status'] != 'complete' or len(ledger['exact_points']) != points:
        raise ValueError('Incomplete curve: '+str(folder))
    expected = dict(electrothermal_linear_solver='umfpack', reuse_linear_analysis=enabled,
                    reuse_sparselu_symbolic=True)
    if ledger['linear_policy'] != expected:
        raise ValueError('Effective default policy mismatch')
    totals = {}
    for section, entries in (('initialization', ledger['initialization_runs']), ('drain', ledger['runs'])):
        row = dict(updates=0, attempts=len(entries), rejected_attempts=0, trials=0, counters={})
        for entry in entries:
            path = entry['result'] if section == 'initialization' else Path(entry['directory'])/'output.json'
            result = read_bound_record(path)
            perf = result['performance']
            if perf['linear_solver'] != 'umfpack' or perf['ialmob_screening_method'] != 'halley':
                raise ValueError('Unexpected effective solver/root method')
            if perf['linear_analysis_reuse_enabled'] != enabled:
                raise ValueError('Unexpected effective reuse policy')
            row['updates'] += result['newton_updates']
            row['trials'] += sum(h.get('line_search_trials', 0) for h in result.get('history', []))
            row['rejected_attempts'] += not entry['gate']['pass_gate']
            for key in COUNTERS:
                row['counters'][key] = row['counters'].get(key, 0)+perf.get(key, 0)
            if section == 'drain' and entry['gate']['pass_gate'] and not state_gate(result, entry['bias_V'])['pass_gate']:
                raise ValueError('Accepted continuation failed original gate')
        totals[section] = row
    for point in ledger['exact_points']:
        if not state_gate(read_bound_record(point['result']), point['bias_V'])['pass_gate']:
            raise ValueError('Exact point failed original gate')
    return ledger, totals


def trajectory(ledger):
    return [(r['bias_V'], r['gate']['pass_gate']) for r in ledger['runs']]


def main():
    if sys.platform != 'linux':
        raise RuntimeError('This closeout driver uses Linux process CPU accounting')
    import resource
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('runner', 'prior', 'native-local', 'output'):
        p.add_argument('--'+key, type=Path, required=True)
    args = p.parse_args()
    runner, prior, native, out = (getattr(args, key).resolve() for key in ('runner', 'prior', 'native_local', 'output'))
    root = Path(__file__).resolve().parents[1]
    out.mkdir(parents=True, exist_ok=False)
    os.environ.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', OMP_DYNAMIC='FALSE', VELA_LINEAR_THREADS='1')
    os.environ.pop('VELA_BLAS_THREADS', None)
    os.environ.pop('VELA_LINEAR_SOLVER', None)  # Exercise compiled default.
    threads = subprocess.check_output([sys.executable, '-c',
        'import ctypes; print(ctypes.CDLL("libopenblas.so.0").openblas_get_num_threads())'], text=True).strip()
    if threads != '1':
        raise ValueError('BLAS thread configuration not one')
    status = dict(status='running', stage='freeze', pid=os.getpid(), started=time.time(),
                  runs=[], comparisons=[], native_checks=[], checkpoints=[], source_commit='ff9fccb1',
                  scope='D0, Vg4/8, neutral 300 K initialization, UMFPACK 1/1, unchanged gates',
                  blas_threads=int(threads))
    publish = lambda: save(out/'summary.json', status)
    frozen = {}
    def freeze(path):
        path = Path(path).resolve()
        if not path.is_file():
            raise ValueError('Missing frozen dependency: '+str(path))
        frozen[str(path)] = digest(path)
    def freeze_dependencies(value):
        if isinstance(value, dict):
            for child in value.values(): freeze_dependencies(child)
        elif isinstance(value, list):
            for child in value: freeze_dependencies(child)
        elif isinstance(value, str) and value.startswith('/') and Path(value).is_file():
            freeze(value)
    for directory in ('src', 'include', 'scripts', 'tests'):
        for path in (root/directory).rglob('*'):
            if path.is_file() and '__pycache__' not in path.parts: freeze(path)
    for path in (runner, runner.parent/'silicon_thermal_probe', root/'CMakeLists.txt', runner.parent/'CMakeCache.txt'):
        freeze(path)
    contract = root/'reference_tcad/templates_ldmos_sentaurus2022/thermal/d0_thermal_acceptance.json'
    local_contract = contract.with_name('d0_joint_acceptance_20260914_v2.json')
    freeze(contract); freeze(local_contract)
    for base in (prior/'evidence/native_20260925', native):
        for path in base.rglob('*'):
            if path.is_file(): freeze(path)
    for gate in (4, 8):
        for name in ('input.json', 'deck.json'):
            source = prior/f'evidence/full_20260925/vg{gate}_halley'/name
            freeze(source)
            freeze_dependencies(read(source))
    save(out/'freeze.json', frozen)
    status['frozen_manifest_sha256'] = digest(out/'freeze.json')
    def verify():
        if any(digest(Path(path)) != value for path, value in frozen.items()):
            raise ValueError('Frozen source, executable, reference or input changed')
    def execute(name, command, expected=0):
        if shutil.disk_usage(out).free < 700*1024**2:
            raise RuntimeError('Less than 700 MiB free; preserve evidence and stop')
        status['current'] = name
        usage = resource.getrusage(resource.RUSAGE_CHILDREN)
        before = usage.ru_utime+usage.ru_stime
        start = time.perf_counter()
        load = os.getloadavg()
        with (out/(name+'.log')).open('x') as log:
            child = subprocess.Popen(list(map(str, command)), cwd=root, stdout=log, stderr=subprocess.STDOUT)
            status['child_pid'] = child.pid; publish()
            rc = child.wait()
        usage = resource.getrusage(resource.RUSAGE_CHILDREN)
        row = dict(name=name, wall_seconds=time.perf_counter()-start,
                   cpu_seconds=usage.ru_utime+usage.ru_stime-before, returncode=rc, load_before=load)
        status.pop('child_pid', None); status['runs'].append(row); publish()
        if rc != expected: raise RuntimeError(f'{name}: exit {rc}, expected {expected}')
        return row
    def prepare(name, gate, enabled, points):
        folder = out/name; folder.mkdir()
        original = prior/f'evidence/full_20260925/vg{gate}_halley'
        save(folder/'input.json', policy_input(read(original/'input.json'), enabled))
        deck = policy_input(read(original/'deck.json'), enabled)
        deck.update(input_file=str(folder/'input.json'), output_directory=str(folder/'results'), resume=False, pause_after_attempts=0)
        deck['sweep']['bias_points_V'] = deck['sweep']['bias_points_V'][:points]
        if deck['initialization']['mode'] != 'neutral_300K': raise ValueError('Neutral initialization required')
        save(folder/'deck.json', deck)
        return folder, deck
    def compare(left, right, scope):
        a, b = (read(f/'results/ledger.json') for f in (left, right))
        if len(a['exact_points']) != len(b['exact_points']): raise ValueError('Point count mismatch')
        for pa, pb in zip(a['exact_points'], b['exact_points']):
            if pa['bias_V'] != pb['bias_V']: raise ValueError('Bias mismatch')
            status['comparisons'].append(dict(scope=scope, bias_V=pa['bias_V'],
                errors=compare_states(read_bound_record(pa['result']), read_bound_record(pb['result']))))
        status.setdefault('trajectories', []).append(dict(scope=scope, equal=trajectory(a)==trajectory(b)))
        publish()
    def score(folders, name):
        result = analyze(prior/'evidence/native_20260925', {g:folders[g]/'results' for g in (4,8)}, read(contract), 15)
        save(out/(name+'.json'), result)
        status['native_checks'].append(dict(name=name, status=result['status'])); publish()
        if result['status'] != 'pass': raise ValueError('Original native electrical/thermal gates failed')
    def curves(stage, points, repetition=0):
        status['stage'] = stage; publish(); verify()
        folders = {}
        for gate in (4,8):
            for enabled in ((False,True) if (repetition+gate//4)%2 else (True,False)):
                name = f'{stage}_vg{gate}_{int(enabled)}'
                folder, deck = prepare(name, gate, enabled, points)
                row = execute(name, [runner, '--config', folder/'deck.json'])
                ledger, totals = inspect(folder, points, enabled)
                row.update(gate=gate, enabled=enabled, points=points, totals=totals, stage=stage)
                folders[gate,enabled] = folder; publish()
                if stage.startswith('repeat'):
                    compare(full[gate,enabled], folder, name+'_vs_full')
            compare(folders[gate,False], folders[gate,True], stage+f'_vg{gate}_toggle')
        if points == 31:
            for enabled in (False,True): score({g:folders[g,enabled] for g in (4,8)}, stage+f'_joint_{int(enabled)}')
        return folders
    publish()
    try:
        curves('short', 8)
        full = curves('full', 31)
        status['stage'] = 'native_local_fields'; publish()
        field_curves = out/'field_curves'; field_curves.mkdir()
        for gate in (4,8): (field_curves/f'vg{gate}').symlink_to(full[gate,True]/'results', target_is_directory=True)
        execute('local_fields', [sys.executable, root/'scripts/audit_templates_ldmos_joint_local_fields.py',
            '--curves',field_curves,'--prefix','vg','--probe',runner.parent/'silicon_thermal_probe',
            '--native',native,'--native-raw',native/'raw','--contract',local_contract,'--output',out/'local_fields'])
        local = read(out/'local_fields/summary.json')
        if local['status'] != 'pass' or len(local['points']) != 62: raise ValueError('Native local-field qualification failed')
        status['native_checks'].append(dict(name='local_fields_62',status='pass')); publish()
        # Checkpoints at known low/mid/high accepted attempts, with a new process
        # at every restart. Do not duplicate or edit previously accepted records.
        status['stage'] = 'recovery'; publish()
        recovered = {}
        for gate in (4,8):
            base = read(full[gate,True]/'results/ledger.json')
            stops = [next(i+1 for i,r in enumerate(base['runs']) if r['gate']['pass_gate'] and r['bias_V']>=bias) for bias in (1.333333333333,20.,36.)]
            folder, deck = prepare(f'recovery_vg{gate}',gate,True,31)
            previous = 0
            for part, stop in enumerate(stops+[0]):
                deck.update(resume=part>0, pause_after_attempts=(stop-previous if stop else 0))
                path=folder/f'deck_{part}.json';save(path,deck)
                row=execute(f'recovery_vg{gate}_part{part}',[runner,'--config',path],expected=1 if stop else 0)
                ledger=read(folder/'results/ledger.json')
                if stop:
                    if ledger['status']!='stopped_at_checkpoint' or len(ledger['runs'])!=stop:
                        raise ValueError('Unexpected checkpoint boundary')
                    status['checkpoints'].append(dict(gate=gate,bias_V=ledger['accepted_bias_V'],attempts=stop,part=part))
                if part:
                    first=read_bound_record(Path(ledger['runs'][previous]['directory'])/'output.json')['performance']
                    if first['linear_object_reused'] or first['symbolic_analyses'] < 1:
                        raise ValueError('Restart must rebuild fresh linear analysis')
                previous=stop;publish()
            inspect(folder,31,True);compare(full[gate,True],folder,f'recovery_vg{gate}')
            recovered[gate]=folder
        score(recovered,'recovery_joint')
        for repetition in range(3): curves(f'repeat{repetition}',31,repetition)
        verify()
        status.update(status='completed',stage='qualified',finished=time.time());publish()
    except BaseException as error:
        status.update(status='failed',error=repr(error),traceback=traceback.format_exc(),finished=time.time());publish();raise


if __name__ == '__main__': main()
