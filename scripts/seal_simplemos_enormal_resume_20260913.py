"""Close the resumed cohort without hiding administrative interruptions."""
from datetime import datetime
from pathlib import Path

import resume_simplemos_enormal_curves_20260913 as resume

e, a, d, O = resume.e, resume.a, resume.d, resume.OUT


def main():
    evidence = [O / name for name in ('preflight_evidence.json', 'completed_evidence.json',
                                    'continuation_resume_evidence.json', 'native_resume_evidence.json')]
    evidence += [e.O / 'completion_evidence.json', e.O / 'fields/visual_evidence.json']
    for path in evidence:
        a.verify(path)
    final = a.read(e.O / 'final_summary.json')
    assert final['all_qualified'] and final['comparison_qualified'] == 204
    assert a.read(e.O / 'fields/summary.json')['all_reconstruction_qualified']
    before = a.rows(O / 'completed_before_resume.csv')
    after = [row for arm in ('continuation', 'native') for row in a.rows(e.O / (arm + '_attempts.csv'))]
    assert len(before) == 306
    old_keys = {(r['case'], int(r['index']), r['arm']) for r in before}
    new_keys = {(r['case'], int(r['index']), r['arm']) for r in after} - old_keys
    assert len(new_keys) == 102
    summary = dict(time=datetime.now().astimezone().isoformat(),
                   completed_before_resume=306, remaining_targets_completed=102,
                   all_targets=408, qualified_bias_pairs=204,
                   completed_attempts=len(after), numerical_failed_attempts=final['failed_attempts'],
                   numerical_reloads=final['reloads'],
                   administrative_interrupted_attempts=len(a.read(O / 'interrupted.json')),
                   previous_306_states_byte_identical=True, interrupted_archives_verified=True,
                   all_qualified=True, acceptance_changed=False, production_changed_during_resume=False,
                   scope='Enormal four complete 0..1 V curves and native field/SRH comparison. Four interrupted attempts are separate from completed numerical attempts. HFS and original full matrix not executed.')
    a.write(O / 'completion_summary.json', summary)
    d.matrix.freeze(O / 'completion_evidence.json',
                    evidence + [Path(__file__).resolve(), O / 'completion_summary.json'])
    print(summary, flush=True)


if __name__ == '__main__':
    main()
