"""Preserve interrupted work and resume the unchanged Enormal curve driver."""
import argparse
import shutil
from datetime import datetime
from pathlib import Path

import simplemos_enormal_curves_20260912 as e

a, d = e.a, e.d
OUT = e.O / 'resume_20260913'
LOCAL = e.L / 'resume_20260913'


def prepare():
    e.configure()
    for name in ('supplemental_source_evidence.json', 'control_evidence.json', 'export_evidence.json'):
        a.verify(e.O / name)
    assert a.read(e.O / 'control_summary.json')['all_qualified']
    assert not (OUT / 'preflight_evidence.json').exists()
    completed, interrupted = [], []
    preserved = []
    for root in sorted((e.L / 'vela').rglob('attempt_*')):
        if not root.is_dir():
            continue
        assert root.resolve().is_relative_to((e.L / 'vela').resolve())
        a.verify(root / 'input_freeze.json')
        if (root / 'result.json').exists():
            result = a.read(root / 'result.json')
            assert result['qualified'] and (root / 'config.status.json').exists()
            completed.append(result)
            preserved += [p for p in root.iterdir() if p.is_file()]
        else:
            assert not (root / 'config.status.json').exists(), 'Interrupted postprocessing needs separate handling'
            destination = LOCAL / 'interrupted' / root.relative_to(e.L / 'vela')
            shutil.copytree(root, destination)
            hashes = {p.name: a.sha(p) for p in root.iterdir() if p.is_file()}
            assert hashes == {p.name: a.sha(p) for p in destination.iterdir() if p.is_file()}
            interrupted.append(dict(original=str(root), archive=str(destination), files=hashes,
                                    classification='Administrative interruption: no final solver status or result; not a qualified state or a numerical failure.'))
    assert len(completed) == 306 and len(interrupted) == 4
    for name in ('continuation.log', 'native.log', 'active_run.json',
                 'postprocessing_status.json', 'postprocessing_coordinator.log'):
        source = e.L / name
        if source.exists():
            destination = LOCAL / 'previous_logs' / name
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, destination)
            assert a.sha(source) == a.sha(destination)
    a.write(OUT / 'interrupted.json', interrupted)
    a.write_csv(OUT / 'completed_before_resume.csv', completed)
    summary = dict(time=datetime.now().astimezone().isoformat(), completed_states=306,
                   remaining_states=102, interrupted_targets=4,
                   runner=str(e.q.RUNNER), runner_sha256=a.sha(e.q.RUNNER),
                   source_changed=False, acceptance_changed=False,
                   method='Original driver and cached completed solver/probe statuses. Repeat four archived unfinished targets from their original input seed; preserve the original one-reload rule for numerical failures.')
    a.write(OUT / 'summary.json', summary)
    d.matrix.freeze(OUT / 'completed_evidence.json', preserved + [OUT / 'completed_before_resume.csv'])
    d.matrix.freeze(OUT / 'preflight_evidence.json',
                    [Path(__file__).resolve(), OUT / 'completed_evidence.json',
                     OUT / 'interrupted.json', OUT / 'summary.json',
                     e.O / 'supplemental_source_evidence.json', e.O / 'control_evidence.json'] +
                    [p for p in LOCAL.rglob('*') if p.is_file()])
    print(summary, flush=True)


def run(arm):
    e.configure()
    a.verify(OUT / 'preflight_evidence.json')
    a.verify(OUT / 'completed_evidence.json')
    a.verify(e.O / 'supplemental_source_evidence.json')
    a.verify(e.O / 'control_evidence.json')
    assert a.read(e.O / 'control_summary.json')['all_qualified']
    e.v.run(arm)
    a.verify(OUT / 'completed_evidence.json')
    d.matrix.freeze(OUT / (arm + '_resume_evidence.json'),
                    [OUT / 'preflight_evidence.json', OUT / 'completed_evidence.json',
                     e.O / ('vela_' + arm + '_evidence.json')])
    print('Finished arm; all 306 previously completed states remain byte-identical:', arm, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=('prepare', 'continuation', 'native'))
    args = parser.parse_args()
    if args.action == 'prepare':
        prepare()
    else:
        run(args.action)
