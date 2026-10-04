"""Seal both bounded validations without replacing the earlier field qualification."""
from pathlib import Path
import re
import validate_simplemos_weak_columns_20260906 as w
import validate_simplemos_minority_invariance_20260906 as m

a=m.a;d=m.d


def main():
    manifests=[m.b.OUT/'followup_evidence.json',w.OUT/'freeze.json',w.OUT/'result.json',
        m.OUT/'freeze.json',m.OUT/'screen_freeze.json',m.OUT/'result.json']
    for path in manifests:a.verify(path)
    files=set(manifests)
    for root in (w.OUT,m.OUT,w.LOCAL,m.LOCAL):
        files.update(p for p in root.rglob('*') if p.is_file()
                     and p.suffix in ('.csv','.json','.log','.txt','.cpp')
                     and p.name!='validation_evidence.json')
    files.add(w.RUNNER)
    scripts=('validate_simplemos_weak_columns_20260906.py','analyze_simplemos_weak_columns_20260906.py',
        'validate_simplemos_minority_invariance_20260906.py','screen_simplemos_minority_invariance_20260906.py',
        'analyze_simplemos_minority_invariance_20260906.py')
    files.update(d.REPO/'scripts'/name for name in scripts)
    files.update([Path(__file__).resolve(),d.REPO/'tests/regression/test_simplemos_weak_and_minority_validation.py'])
    report=d.REPO/'docs/validation/simplemos_weak_columns_and_minority_validation_2026-09-06.md';files.add(report)
    for target in re.findall(r'\]\(([^)]+)\)',report.read_text(encoding='utf-8')):
        if '://' in target:continue
        path=Path(target.split('#')[0]);path=path if path.is_absolute() else report.parent/path
        if path.name!='validation_evidence.json':assert path.exists(),path
    weak=a.read(w.OUT/'result.json');minority=a.read(m.OUT/'result.json')
    assert weak['targets']==weak['reference_qualified']==33
    cases=a.rows(m.OUT/'cases.csv');runs=a.rows(m.OUT/'qualification.csv')
    assert len(cases)==8
    assert len(runs)+minority['not_completed_states']==32
    for row in cases:
        if row['invariance_qualified']=='True':
            assert int(row['initializations'])==int(row['qualified_states'])==4
            assert row['numerical_invariance']=='True'
    assert all(int(r['all_row_checked_rows'])==int(r['expected_active_rows']) for r in runs)
    a.write(m.OUT/'validation_evidence.json',dict(schema_version=1,date='2026-09-06',
        status='completed_bounded_validation_with_failures_retained',
        input_hashes={a.rel(p):a.sha(p) for p in sorted(files)},verified_manifests=len(manifests),
        weak_targets=33,weak_reference_qualified=weak['reference_qualified'],
        weak_original_step_passed=weak['legacy_step_passed'],weak_fine_step_passed=weak['fine_step_passed'],
        minority_completed_states=minority['completed_states'],minority_all_row_qualified=minority['all_row_qualified'],
        minority_invariance_qualified_cases=minority['invariance_qualified_cases'],
        interrupted_attempts=len(a.read(m.OUT/'interrupted_processes.json')['processes']),
        checks={'isolated_cpp_builds':1,'python_regression_tests':7,'python_regression_passed':7,
            'regression_command':'D:/msys64/ucrt64/bin/python.exe tests/regression/test_simplemos_weak_and_minority_validation.py',
            'python_syntax':'Passed for the new validation scripts','report_links':'Checked',
            'historical_inputs':'Prior followup evidence and frozen inputs verified without alteration'},
        limitations=['Local column audit is not a proof over every device column and branch.',
            'Underlying edge flux and mobility evaluations remain double precision.',
            'Failed all-row states cannot qualify initialization invariance or establish physical nonuniqueness.',
            'High NWell current discrepancy is not resolved by this validation.'],
        production_solver_changed=False,old_gate_changed=False,m82_released=False,m83_released=False))
    a.verify(m.OUT/'validation_evidence.json')
    print('Sealed',len(files),'files;',len(manifests),'input manifests verified',flush=True)


if __name__=='__main__':main()
