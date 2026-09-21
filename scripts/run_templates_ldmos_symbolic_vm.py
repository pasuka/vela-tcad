"""Serial same-VM R10/explicit-HFS/native full repeat matrix; original gates."""
import argparse
import calendar
import copy
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import time

from run_templates_ldmos_contact_sweep import read, sha
from run_templates_ldmos_electrothermal_curve import state_gate
from analyze_templates_ldmos_predictor_study import state_delta


def audit(directory, reference):
    ledger, old = read(directory/'results/ledger.json'), read(reference/'results/ledger.json')
    assert ledger['status'] == 'complete' and len(ledger['exact_points']) == 31
    maximum = [0.]*4
    for p, q in zip(ledger['exact_points'], old['exact_points']):
        assert p['bias_V'] == q['bias_V']
        value, frozen = read(p['result']), read(q['result'])
        assert state_gate(value, p['bias_V'])['pass_gate']
        maximum = list(map(max, maximum, state_delta(value, frozen)))
    assert all(d <= limit for d, limit in zip(maximum, (1e-8, 1e-8, 1e-8, 1e-7))), maximum
    assert len(ledger['initialization_runs']) == len(old['initialization_runs'])
    init_max = [0.]*4
    for a, b in zip(ledger['initialization_runs'], old['initialization_runs']):
        init_max = list(map(max, init_max, state_delta(read(a['result']), read(b['result']))))
    assert all(d <= limit for d, limit in zip(init_max, (1e-8, 1e-8, 1e-8, 1e-7))), init_max
    costs = dict(drain_updates=0, initialization_updates=0, trials=0, failed_attempts=0, performance={})
    signatures = []
    for initialization in (True, False):
        for row in ledger['initialization_runs' if initialization else 'runs']:
            value = read(row['result'] if initialization else Path(row['directory'])/'output.json')
            if not initialization:
                assert state_gate(value, row['bias_V'])['pass_gate'] == row['gate']['pass_gate']
                costs['failed_attempts'] += not row['gate']['pass_gate']
            costs['initialization_updates' if initialization else 'drain_updates'] += value['newton_updates']
            costs['trials'] += sum(h['line_search_trials'] for h in value['history'])
            assert value['performance']['linear_solver'] == 'umfpack'
            for key, v in value['performance'].items():
                if isinstance(v, (int, float)):
                    costs['performance'][key] = costs['performance'].get(key, 0)+v
            signatures.append(sha(Path(row['result']) if initialization else Path(row['directory'])/'output.json'))
    return dict(pass_gate=True, max_state_delta=maximum, init_max_state_delta=init_max,
                costs=costs, ledger_sha256=sha(directory/'results/ledger.json'), output_sha256=signatures)


def carry_completed(previous, runner_hash, profile_hash, repeats, envroot, profile, native):
    """Validate immutable completed runs; interrupted costs remain separate."""
    prior = read(previous/'summary.json')
    assert prior['status'] in ('paused_at_deadline', 'paused_before_next_run', 'complete')
    assert prior['runner_sha256'] == runner_hash == sha(previous/'vela_example_runner')
    assert prior['profile_sha256'] == profile_hash and prior['repeats'] == repeats
    completed, interrupted, seen = [], [], set()
    for row in prior['runs']:
        key = (row['repeat'], row['gate_V'], row['variant'])
        assert key not in seen and 0 <= key[0] < repeats and key[1] in (4, 8)
        assert key[2] in ('native', 'baseline', 'high')
        seen.add(key)
        if row.get('interrupted_at_deadline') or row['exit_code'] != 0:
            interrupted.append(copy.deepcopy(row))
            continue
        directory = Path(row['directory'])
        if row['variant'] == 'native':
            for name in ('IdVd.cmd', 'n1_fps.tdr', 'sdevice.par', 'Siliconc100.par'):
                assert sha(directory/name) == native['cases_sha256']['native_vg%d/%s' % (row['gate_V'], name)]
            assert len(list(directory.glob('field_vg%d_*_des.tdr' % row['gate_V']))) == 31
        else:
            assert row.get('qualification', {}).get('pass_gate')
            for name in ('input', 'deck'):
                assert sha(directory/(name+'.json')) == row[name+'_sha256']
            case = next(c for c in profile['cases'] if c['gate_V'] == row['gate_V'])
            assert audit(directory, envroot/Path(case['deck']).parent) == row['qualification']
        completed.append(copy.deepcopy(row))
    return completed, copy.deepcopy(prior.get('prior_interrupted_runs', []))+interrupted


def main():
    import resource  # Linux process accounting; helpers remain testable on Windows.
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('environment', 'profile', 'runner', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--repeats', type=int, default=2)
    parser.add_argument('--continue-from', type=Path, help='Validate/carry completed runs into a NEW output directory; rerun interruptions cold')
    parser.add_argument('--deadline-utc', help='Stop this matrix by YYYY-MM-DDTHH:MM:SSZ; never resume automatically')
    args = parser.parse_args()
    deadline = (calendar.timegm(time.strptime(args.deadline_utc, '%Y-%m-%dT%H:%M:%SZ'))
                if args.deadline_utc else None)
    if args.repeats < 2:
        parser.error('At least two rounds required')
    envroot, out = args.environment.resolve(), args.output.resolve()
    profile = read(args.profile)
    for name, digest in profile['files_sha256'].items():
        assert sha(envroot/name) == digest, name
    out.mkdir(parents=True, exist_ok=False)
    binary = out/'vela_example_runner'
    shutil.copy2(str(args.runner), str(binary))
    native = read(envroot/'benchmark_r7.json')
    report = dict(status='running', runner_sha256=sha(binary), profile_sha256=sha(args.profile),
                  repeats=args.repeats, runs=[], scope=__doc__, deadline_utc=args.deadline_utc)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')

    def save():
        temporary = out/'summary.tmp'
        temporary.write_text(json.dumps(report, indent=2))
        os.replace(str(temporary), str(out/'summary.json'))

    save()
    try:
        if args.continue_from:
            previous = args.continue_from.resolve()
            completed, interrupted = carry_completed(previous, report['runner_sha256'],
                report['profile_sha256'], args.repeats, envroot, profile, native)
            report.update(runs=completed, prior_interrupted_runs=interrupted,
                          continuation_summary=str(previous/'summary.json'),
                          continuation_summary_sha256=sha(previous/'summary.json'))
            save()
            print(json.dumps(dict(carried_completed=len(completed), retained_interruptions=len(interrupted))), flush=True)
        done = {(r['repeat'], r['gate_V'], r['variant']) for r in report['runs']}
        for repeat in range(args.repeats):
            for gate in (4, 8):
                case = next(c for c in profile['cases'] if c['gate_V'] == gate)
                order = ('native', 'baseline', 'high') if (repeat+(gate == 8)) % 2 == 0 else ('high', 'baseline', 'native')
                for variant in order:
                    if (repeat, gate, variant) in done:
                        continue
                    if deadline is not None and time.time() >= deadline-120:
                        report.update(status='paused_before_next_run', next_run=dict(repeat=repeat, gate_V=gate, variant=variant))
                        save()
                        return
                    directory = out/('r%d_vg%d_%s' % (repeat, gate, variant))
                    directory.mkdir()
                    row = dict(repeat=repeat, gate_V=gate, variant=variant, directory=str(directory),
                               start_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
                               load_before=os.getloadavg(), cpu_stat_before=Path('/proc/stat').read_text().splitlines()[0])
                    if variant == 'native':
                        for name in ('IdVd.cmd', 'n1_fps.tdr', 'sdevice.par', 'Siliconc100.par'):
                            source = envroot/'cases_r7'/('native_vg%d' % gate)/name
                            assert sha(source) == native['cases_sha256']['native_vg%d/%s' % (gate, name)]
                            shutil.copy2(str(source), str(directory/name))
                        argv = ['/atctools/Synopsys/tcad/T-2022.03/bin/sdevice', 'IdVd.cmd']
                    else:
                        cfg, deck = copy.deepcopy(read(envroot/case['input'])), copy.deepcopy(read(envroot/case['deck']))
                        assert deck['reuse_static_preparation'] and cfg['reuse_ialmob_thermal_high_field']
                        assert deck['initialization']['mode'] == 'neutral_300K'
                        assert len(deck['sweep']['bias_points_V']) == 31
                        assert not cfg.get('diagnostic_ialmob_kernel_timing', False)
                        cfg['diagnostic_ialmob_explicit_high_field'] = variant == 'high'
                        cfg['diagnostic_ialmob_generated_low_field'] = False
                        deck.update(input_file=str(directory/'input.json'), output_directory=str(directory/'results'))
                        (directory/'input.json').write_text(json.dumps(cfg))
                        (directory/'deck.json').write_text(json.dumps(deck, indent=2))
                        row.update(input_sha256=sha(directory/'input.json'), deck_sha256=sha(directory/'deck.json'))
                        argv = [str(binary), '--config', str(directory/'deck.json')]
                    before = resource.getrusage(resource.RUSAGE_CHILDREN)
                    started = time.perf_counter()
                    with (directory/'run.log').open('x') as log:
                        proc = subprocess.Popen(argv, cwd=str(directory), env=env, stdout=log, stderr=subprocess.STDOUT,
                                                start_new_session=True)
                        report['active_run'] = dict(row, pid=proc.pid)
                        save()
                        paused = False
                        terminated = False
                        try:
                            while True:
                                now = time.time()
                                if deadline is not None and now >= deadline-60:
                                    paused = True
                                    if variant != 'native':
                                        result_root = directory/'results'
                                        result_root.mkdir(exist_ok=True)
                                        (result_root/'STOP').touch()
                                    if now >= deadline-5 and proc.poll() is None and not terminated:
                                        os.killpg(proc.pid, signal.SIGTERM)
                                        terminated = True
                                    if now >= deadline and proc.poll() is None:
                                        os.killpg(proc.pid, signal.SIGKILL)
                                try:
                                    code = proc.wait(timeout=1 if paused else 5)
                                    break
                                except subprocess.TimeoutExpired:
                                    status = dict(run=directory.name, elapsed=time.perf_counter()-started,
                                                  pause_requested=paused)
                                    path = directory/'results/ledger.json'
                                    if path.exists():
                                        ledger = read(path)
                                        status.update(bias_V=ledger['accepted_bias_V'], exact_points=len(ledger['exact_points']))
                                    print(json.dumps(status), flush=True)
                        except BaseException:
                            if proc.poll() is None:
                                os.killpg(proc.pid, signal.SIGTERM)
                                try:
                                    proc.wait(timeout=3)
                                except subprocess.TimeoutExpired:
                                    os.killpg(proc.pid, signal.SIGKILL)
                                    proc.wait()
                            raise
                    wall = time.perf_counter()-started
                    after = resource.getrusage(resource.RUSAGE_CHILDREN)
                    row.update(exit_code=code, external_wall_seconds=wall,
                               child_cpu_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,
                               load_after=os.getloadavg(), cpu_stat_after=Path('/proc/stat').read_text().splitlines()[0])
                    report['runs'].append(row)
                    report.pop('active_run', None)
                    save()  # Failures and audit errors retain their entire process cost.
                    if paused:
                        row['interrupted_at_deadline'] = True
                        report['status'] = 'paused_at_deadline'
                        save()
                        return
                    assert code == 0, row
                    if variant != 'native':
                        row['qualification'] = audit(directory, envroot/Path(case['deck']).parent)
                    save()
                    print(json.dumps(row), flush=True)
        assert sha(binary) == report['runner_sha256']
        report['status'] = 'complete'
        save()
    except BaseException as error:
        report.update(status='failed', error=repr(error))
        save()
        raise


if __name__ == '__main__':
    main()
