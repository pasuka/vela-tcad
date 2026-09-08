"""Seal the bounded refinement A/B without changing historical failures or gates."""
from pathlib import Path
import ast
import re
import validate_simplemos_linear_refinement_20260906 as m

a=m.a;d=m.d


def main():
    manifests=[m.v.OUT/'validation_evidence.json',m.OUT/'freeze.json',m.OUT/'result.json',m.OUT/'detail_result.json']
    for p in manifests:a.verify(p)
    contract=a.read(m.OUT/'contract.json');result=a.read(m.OUT/'result.json');rows=a.rows(m.OUT/'comparison.csv')
    assert result['completed_states']==len(rows)==7
    for j in contract['disabled_sentinels']:
        p=Path(j['config']).parent;src=Path(j['source']);s=a.read(p/'config.status.json');old=a.read(src/'config.status.json')
        assert a.sha(p/'state.csv')==a.sha(src/'state.csv')
        assert all(s[k]==old[k] for k in ('iterations','converged','exit_code','failure_reason'))
    for r in rows:
        if r['refined_qualified']=='True':
            assert r['refined_converged']=='True' and int(r['refined_violations'])==int(r['zero_scale_rows'])==0
            assert r['all_row_global_satisfied']=='True' and float(r['refined_kcl_over_Id'])<=1e-8
    tests=(m.LOCAL/'test_result.log').read_text();assert 'All tests passed (17 assertions in 4 test cases)' in tests
    files=set(manifests)
    for root in (m.OUT,m.LOCAL):
        files.update(p for p in root.rglob('*') if p.is_file() and p.suffix in ('.csv','.json','.txt','.log','.cpp','.exe') and p.name!='validation_evidence.json')
    scripts=[d.REPO/'scripts'/n for n in ('validate_simplemos_linear_refinement_20260906.py','analyze_simplemos_linear_refinement_20260906.py','detail_simplemos_refinement_failures_20260906.py','seal_simplemos_linear_refinement_20260906.py')]
    for p in scripts:ast.parse(p.read_text(encoding='utf-8-sig'),filename=str(p))
    files.update(scripts);files.update((m.HEADER,m.TEST))
    report=d.REPO/'docs/validation/simplemos_linear_refinement_self_consistent_2026-09-06.md';files.add(report)
    text=report.read_text(encoding='utf-8');assert '<!--' not in text
    for target in re.findall(r'\]\(([^)]+)\)',text):
        if '://' in target:continue
        p=Path(target.split('#')[0]);p=p if p.is_absolute() else report.parent/p
        if p.name!='validation_evidence.json':assert p.exists(),p
    a.write(m.OUT/'validation_evidence.json',dict(date='2026-09-06',schema_version=1,status='completed_bounded_self_consistent_refinement_ab',
        input_hashes={a.rel(p):a.sha(p) for p in sorted(files)},verified_manifests=len(manifests),
        disabled_sentinel_replays=2,refined_states=7,baseline_qualified=result['baseline_qualified'],refined_qualified=result['refined_qualified'],
        newly_qualified=result['newly_qualified'],lost_qualification=result['lost_qualification'],
        linear_steps=result['linear_steps'],linear_monitor_passed_steps=result['linear_monitor_passed_steps'],qualified_four_seed_cases=result['qualified_four_seed_cases'],
        checks=dict(isolated_solver_build='passed',catch2_tests=4,catch2_assertions=17,catch2_result='passed',
            python_syntax_files=len(scripts),report_links='passed',historical_manifest='verified',full_ctest='not_run_isolated_diagnostic_only'),
        production_solver_changed=False,physical_model_changed=False,nonlinear_gate_changed=False,line_search_changed=False,remote_runs=0,
        m82_released=False,m83_released=False,
        limitations=['Only three bias points and seven selected original initial states.',
            'The linear accuracy monitor is not a substitute for nonlinear qualification.',
            'Failed states remain unsuitable for physical field attribution or initialization-invariance certification.',
            'No full Id-Vg revalidation or production algorithm change in this experiment.']))
    a.verify(m.OUT/'validation_evidence.json');print('Sealed',len(files),'files;',len(manifests),'manifests verified',flush=True)


if __name__=='__main__':main()
