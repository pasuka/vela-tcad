"""Seal fixed-state arithmetic and update-localization evidence, retaining failed states."""
from pathlib import Path
import ast
import re
import audit_simplemos_minority_residual_20260906 as v
import check_simplemos_minority_linear_20260906 as l

a=v.a;d=v.d


def main():
    manifests=[v.m.OUT/'validation_evidence.json',v.OUT/'freeze.json',v.OUT/'precision_result.json',
        v.OUT/'trace_result.json',l.OUT/'freeze.json',l.OUT/'result.json']
    for p in manifests:a.verify(p)
    files=set(manifests)
    for root in (v.OUT,v.LOCAL):
        files.update(p for p in root.rglob('*') if p.is_file() and p.suffix in
            ('.csv','.json','.log','.txt','.cpp','.h','.exe') and p.name!='validation_evidence.json')
    scripts=('audit_simplemos_minority_residual_20260906.py','analyze_simplemos_minority_residual_20260906.py',
        'analyze_simplemos_minority_trace_20260906.py','check_simplemos_minority_linear_20260906.py',
        'analyze_simplemos_minority_linear_20260906.py','seal_simplemos_minority_residual_20260906.py')
    python_files=[d.REPO/'scripts'/name for name in scripts]+[d.REPO/'tests/regression/test_simplemos_minority_residual_precision.py']
    for p in python_files:ast.parse(p.read_text(encoding='utf-8-sig'),filename=str(p))
    files.update(python_files)
    report=d.REPO/'docs/validation/simplemos_minority_residual_and_update_localization_2026-09-06.md';files.add(report)
    prose=report.read_text(encoding='utf-8');assert '<!--' not in prose
    for target in re.findall(r'\]\(([^)]+)\)',prose):
        if '://' in target:continue
        p=Path(target.split('#')[0]);p=p if p.is_absolute() else report.parent/p
        if p.name!='validation_evidence.json':assert p.exists(),p
    precision=a.read(v.OUT/'precision_result.json');rows=a.rows(v.OUT/'precision_summary.csv')
    assert precision['states']==len(rows)==20 and precision['rows']==36280
    assert all(r['baseline_violations']==r['kernel_violations']==r['partition_violations'] for r in rows)
    qualified=sum(int(r['baseline_violations'])==0 for r in rows);assert qualified==2
    traces=a.rows(v.OUT/'traces.csv');identity=[r for r in traces if r['mode']=='identity_replay']
    assert len(identity)==6 and all(r['state_identical']=='True' for r in identity)
    single=[r for r in traces if r['mode']=='one_step_saved_failure'];assert len(single)==1 and single[0]['converged']=='False'
    cfg=a.read(Path(single[0]['config']));assert cfg['solver']['max_iter']==cfg['solver']['carrier_row_convergence']['min_newton_max_iter']==1
    trace=a.read(v.OUT/'trace_result.json');assert trace['trials']==70 and trace['all_trial_norms_reconstructed']
    linear=a.read(l.OUT/'result.json');assert linear['states']==3 and linear['groups']==30 and linear['rows']==54420
    assert linear['all_legacy_rows_and_raw_steps_identical'] and linear['nonlinear_solves']==0
    # Hash core files as a direct entry too; prior frozen source manifests remain authoritative.
    files.update(d.REPO/p for p in ('src/equation/CoupledDDAssembler.cpp','src/physics/MobilityModel.cpp','src/solver/NewtonSolver.cpp','src/solver/LinearSolver.cpp','src/numerics/LineSearch.cpp','src/tools/vela_example_runner.cpp'))
    assert all(p.exists() for p in files)
    evidence=dict(schema_version=1,date='2026-09-06',status='completed_localization_with_original_failures_retained',
        input_hashes={a.rel(p):a.sha(p) for p in sorted(files)},verified_manifests=len(manifests),
        fixed_states=20,active_rows=36280,arithmetic_reference_changed_violation_counts=0,all_row_qualified_states=qualified,
        exact_replays=6,single_step_saved_failure=1,reconstructed_last_trials=70,
        static_linear_states=3,static_linear_groups=30,static_linear_active_rows=54420,static_trial_is_dc_qualification=False,
        checks=dict(isolated_cpp_builds_passed=2,initial_static_build_header_failure_retained=True,
            precision_property_tests_run=6,precision_property_tests_passed=6,
            test_command='D:/msys64/ucrt64/bin/python.exe tests/regression/test_simplemos_minority_residual_precision.py',
            python_syntax_files=len(python_files),report_links='passed',prior_seal_and_freezes='verified',
            full_ctest='not_run_no_production_changes'),
        conclusions=['Fixed-state SG/SRH arithmetic changes are too small to remove any original violating rows.',
            'n19 hole node 41 has reproducible weak-row linear solve error; same-matrix refinement restores an effective update.',
            'Current scalar line-search merit is dominated by Poisson in the late failed n19 trajectory.',
            'High-Vd full raw update remains nonlinear after linear defects are corrected; equilibration alone is not uniformly better.'],
        limitations=['Mobility and lifetime coefficients are kept at production double precision.',
            'Linear controls use the existing complete diagnostic Jacobian; they do not independently prove arbitrary Jacobian derivatives.',
            'Static trial ratios use frozen original row scales and do not qualify DC convergence or initialization invariance.',
            'The high NWell Id-Vg discrepancy remains unresolved.'],
        production_solver_changed=False,old_gate_changed=False,remote_runs=0,m82_released=False,m83_released=False)
    a.write(v.OUT/'validation_evidence.json',evidence);a.verify(v.OUT/'validation_evidence.json')
    print('Sealed',len(files),'files; verified',len(manifests),'manifests;',len(python_files),'Python syntax checks',flush=True)


if __name__=='__main__':main()
