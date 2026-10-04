"""Resume the frozen HFS baseline in a separate tree; retain the pause evidence."""
import argparse
import copy
import json
import os
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
import simplemos_hfs_curves_20260914 as e

a, d, v = e.a, e.d, e.v
PAUSE = e.O/'pause_20260914_2213'
LOCAL = e.L/'resume_20260915'
OUT = e.O/'resume_20260915'
ORIGINAL_SOLVE = v.solve_target


def key(row):
    return row['case'], int(row['index']), row['arm']


def load_cache():
    rows = a.read(PAUSE/'completed_states.json')
    assert len(rows) == 64 and all(r['qualified'] and r['attempt'] == 0 for r in rows)
    cache = {key(r): r for r in rows}
    assert len(cache) == 64
    return cache


def check_cached(c, index, arm, seed, row):
    dest = Path(row['dest'])
    assert a.read(dest/'result.json') == row
    expected = copy.deepcopy(a.read(Path(c['template'])))
    for contact in expected['contacts']:
        if contact['name'] == 'gate':
            contact['bias'] = e.n.GRID[index]
    expected.update(state_file=str(seed), output_state_file=str(dest/'state.csv'))
    assert a.read(dest/'config.json') == expected, (c['case'], index, arm)
    a.verify(dest/'input_freeze.json')


def preflight():
    e.configure(); e.gate()
    names = ['supplemental_source_evidence.json', 'control_evidence.json',
             'export_evidence.json', 'sixteen_evidence.json', 'curve_execution_freeze.json',
             'vela_freeze.json', 'pause_20260914_2213/snapshot_evidence.json']
    for name in names:
        a.verify(e.O/name)
    assert not LOCAL.exists() and not OUT.exists()
    for arm in ('continuation', 'native'):
        assert not (e.O/(arm+'_attempts.csv')).exists()
        assert not (e.O/('vela_'+arm+'_evidence.json')).exists()
    cache = load_cache()
    interrupted = {key(r): r for r in a.read(PAUSE/'interrupted_attempts.json')}
    assert len(interrupted) == 4 and not (set(cache) & set(interrupted))
    schedule = []

    def dry_solve(c, index, arm, seed):
        k = c['case'], index, arm
        if k in cache:
            check_cached(c, index, arm, seed, cache[k])
            result = cache[k]
        else:
            result = dict(case=c['case'], index=index, arm=arm, qualified=True,
                          dest=str(LOCAL/'vela'/c['case']/arm/f'vg_{index:03d}'/'attempt_0'))
        if k in interrupted:
            assert Path(interrupted[k]['state_file']) == Path(seed), k
        schedule.append(dict(case=c['case'], index=index, arm=arm, cached=k in cache,
                             previous_interruption=k in interrupted, seed=str(seed), dest=result['dest']))
        return [result]

    # Exercise the actual scheduling functions without solving or creating seeds.
    old_solve, old_seed = v.solve_target, v.native_seed
    try:
        v.solve_target = dry_solve
        v.native_seed = lambda c, i: e.L/'initial'/c['case']/f'vg_{i:03d}'/'state.csv'
        for c in a.read(e.O/'vela_contract.json')['cases']:
            v.continuation(c); v.native_path(c)
    finally:
        v.solve_target, v.native_seed = old_solve, old_seed
    assert len(schedule) == len({key(r) for r in schedule}) == 408
    assert sum(r['cached'] for r in schedule) == 64
    assert sum(r['previous_interruption'] for r in schedule) == 4
    a.write(OUT/'schedule.json', schedule)
    a.write(OUT/'contract.json', dict(cached_states=64, remaining_states=344,
        administrative_interruptions=4, numerical_failures_before_resume=0,
        workers_per_arm=2, arms=['continuation', 'native'], gates=v.GATES,
        original_results_immutable=True, source_and_physics_changed=False,
        output_root=str(LOCAL), runner=str(e.q.RUNNER), runner_sha256=a.sha(e.q.RUNNER),
        recovery='Original maximum one same-bias reload; qualified seeds only.',
        scope='Complete unchanged HFS baseline before SRH volume repair.'))
    d.matrix.freeze(OUT/'preflight_evidence.json', [e.O/name for name in names]+[
        OUT/'contract.json', OUT/'schedule.json', Path(__file__).resolve(),
        Path(__file__).with_name('status_simplemos_hfs_resume_20260915.py'),
        PAUSE/'completed_states.json', PAUSE/'interrupted_attempts.json'])
    LOCAL.mkdir(parents=True)
    print('Preflight passed: 64 immutable cached states, 344 new targets, 4 isolated reruns.', flush=True)


def run(arm):
    e.configure(); e.gate()
    a.verify(OUT/'preflight_evidence.json')
    a.verify(PAUSE/'snapshot_evidence.json')
    cache = load_cache()
    interrupted = {key(r): r for r in a.read(PAUSE/'interrupted_attempts.json')}
    cases = a.read(e.O/'vela_contract.json')['cases']
    assert not (e.O/(arm+'_attempts.csv')).exists()
    assert all(not (LOCAL/'vela'/c['case']/arm).exists() for c in cases)
    marker = LOCAL/(arm+'_active.json')
    assert not marker.exists(), marker

    def mark(status):
        marker.write_text(json.dumps(dict(pid=os.getpid(), arm=arm, status=status,
            recorded_at=datetime.now().astimezone().isoformat()), indent=2), encoding='utf-8')

    def resumed_solve(c, index, selected_arm, seed):
        k = c['case'], index, selected_arm
        if k in cache:
            check_cached(c, index, selected_arm, seed, cache[k])
            return [cache[k]]
        if k in interrupted:
            assert Path(seed) == Path(interrupted[k]['state_file'])
        return ORIGINAL_SOLVE(c, index, selected_arm, seed)

    v.LOCAL = LOCAL
    v.solve_target = resumed_solve
    mark('running')
    try:
        with ThreadPoolExecutor(max_workers=2) as pool:
            rows = [r for group in pool.map(v.continuation if arm == 'continuation' else v.native_path, cases) for r in group]
        assert len({key(r) for r in rows}) == 204
        a.verify(PAUSE/'snapshot_evidence.json')
        v.csv_union(e.O/(arm+'_attempts.csv'), rows)
        files = [OUT/'preflight_evidence.json', PAUSE/'snapshot_evidence.json',
                 e.O/'vela_freeze.json', e.O/(arm+'_attempts.csv')]
        files += [p for c in cases for p in (LOCAL/'vela'/c['case']/arm).rglob('*') if p.is_file()]
        d.matrix.freeze(e.O/('vela_'+arm+'_evidence.json'), files)
        mark('complete')
    except BaseException:
        mark('failed_or_interrupted')
        raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=['preflight', 'continuation', 'native'])
    action = parser.parse_args().action
    preflight() if action == 'preflight' else run(action)
