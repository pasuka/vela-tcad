"""Seal the completed baseline, retaining administrative interruptions separately."""
from pathlib import Path
from datetime import datetime
import resume_simplemos_hfs_curves_20260915 as resume

e, a, d = resume.e, resume.a, resume.d


def main():
    evidence = [resume.OUT/'preflight_evidence.json', resume.OUT/'postprocessing_evidence.json', resume.PAUSE/'snapshot_evidence.json',
                e.O/'vela_continuation_evidence.json', e.O/'vela_native_evidence.json',
                e.O/'completion_evidence.json', e.O/'fields/evidence.json',
                e.O/'srh_volume_evidence.json']
    for path in evidence:
        a.verify(path)
    final = a.read(e.O/'final_summary.json')
    assert final['all_qualified'] and final['comparison_qualified'] == 204
    assert a.read(e.O/'fields/summary.json')['all_reconstruction_qualified']
    before = resume.load_cache()
    after = [r for arm in ('continuation', 'native') for r in a.rows(e.O/(arm+'_attempts.csv'))]
    after_keys = {resume.key(r) for r in after}
    assert len(after_keys) == 408 and len(after_keys-set(before)) == 344
    for k, row in before.items():
        assert a.read(Path(row['dest'])/'result.json') == row
        selected = [r for r in after if resume.key(r) == k]
        assert len(selected) == 1 and selected[0]['dest'] == row['dest']
    for row in after:
        if resume.key(row) not in before:
            assert Path(row['dest']).is_relative_to(resume.LOCAL/'vela')
    native = a.rows(e.O/'comparison.csv')
    worst = max(native, key=lambda r: abs(float(r['continuation_error_percent'])))
    result = dict(recorded_at=datetime.now().astimezone().isoformat(),
        old_completed_states_byte_identical=64, remaining_targets_completed=344,
        all_targets=408, qualified_bias_pairs=204,
        numerical_failed_attempts=final['failed_attempts'], numerical_reloads=final['reloads'],
        administrative_interrupted_attempts=len(a.read(resume.PAUSE/'interrupted_attempts.json')),
        administrative_interruption_files_byte_identical=True,
        all_numerical_and_initialization_gates_qualified=True,
        worst_current_error=dict(case=worst['case'], vg=float(worst['vg']),
            error_percent=float(worst['continuation_error_percent'])),
        current_agreement_within_one_percent=all(abs(float(r['continuation_error_percent']))<=1 for r in native),
        acceptance_changed=False, source_or_physics_changed=False,
        scope='Frozen HFS baseline only. SRH remains all_cell; volume repair and original full matrix remain later work.')
    a.write(resume.OUT/'completion_summary.json', result)
    d.matrix.freeze(resume.OUT/'completion_evidence.json', evidence+[
        Path(__file__).resolve(), resume.OUT/'completion_summary.json'])
    print(result, flush=True)


if __name__ == '__main__':
    main()
