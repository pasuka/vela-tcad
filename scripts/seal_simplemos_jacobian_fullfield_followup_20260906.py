"""Verify frozen inputs and seal this turn's completed numerical investigations."""
from pathlib import Path
import json
import re
import run_simplemos_fullfield_vela_20260906 as b
import audit_simplemos_jacobian_columns_20260906 as columns
import audit_simplemos_surface_derivative_step_20260906 as surface
import analyze_simplemos_channel_poisson_coupling_20260906 as poisson
import refine_simplemos_native_precision_20260906 as precision

a=b.a;d=b.d


def main():
    manifests=[columns.OUT/'freeze.json',columns.OUT/'result.json',surface.OUT/'freeze.json',surface.OUT/'result.json',poisson.OUT/'freeze.json']
    manifests += [b.OUT/(name+'.json') for name in ('native_freeze','export_contract','vela_freeze','recovery_freeze','reclosure_freeze','analysis')]
    manifests += [b.OUT/'visibility/result.json',precision.OUT/'freeze.json',precision.OUT/'export_contract.json',precision.OUT/'precision_result.json',
                  precision.OUT/'comparison/analysis.json',precision.OUT/'comparison/visibility/result.json']
    for path in manifests:a.verify(path)
    files=set(manifests)
    for root in (columns.OUT,poisson.OUT,b.OUT):
        files.update(p for p in root.rglob('*') if p.is_file() and p.name!='followup_evidence.json')
    for root in (columns.LOCAL,poisson.LOCAL):
        files.update(p for p in root.rglob('*') if p.is_file() and p.suffix in ('.csv','.json','.log','.txt','.cpp'))
    for c in a.read(poisson.v.CONTRACT)['cases']:
        files.update([Path(c['mapped']),d.matrix.spatial.precision.LOCAL/c['case']/'state.csv'])
        files.update(d.LOCAL/c['case']/role/'residual.csv' for role in ('strict','mapped'))
    files.update([b.LOCAL/'comparison_nodes.csv',precision.LOCAL/'comparison_nodes.csv',
        d.REPO/'tests/regression/test_simplemos_jacobian_fullfield_followup.py',Path(__file__).resolve()])
    report=d.REPO/'docs/validation/simplemos_jacobian_and_fullfield_followup_2026-09-06.md';files.add(report)
    for target in re.findall(r'\]\(([^)]+)\)',report.read_text(encoding='utf-8')):
        if '://' in target:continue
        path=Path(target.split('#')[0]);path=path if path.is_absolute() else report.parent/path
        if path.name!='followup_evidence.json':assert path.exists(),path
    points=a.rows(precision.OUT/'comparison/points.csv')
    assert len(points)==84 and all(r['qualified']=='True' for r in points)
    assert all(float(r['vela_kcl_over_Id'])<=1e-8 and float(r['native_kcl_over_Id'])<=1e-8 for r in points)
    assert a.read(surface.OUT/'result.json')['stable_failures_remaining']==0
    assert a.read(precision.OUT/'precision_result.json')['kcl_qualified']==42
    # The only solver-core copies altered by this investigation live under build-release.
    for path in (d.REPO/'src/equation/CoupledDDAssembler.cpp',d.REPO/'src/physics/MobilityModel.cpp',d.REPO/'src/solver/NewtonSolver.cpp',d.REPO/'src/tools/vela_example_runner.cpp'):
        files.add(path)
    a.write(b.OUT/'followup_evidence.json',dict(schema_version=1,date='2026-09-06',
        status='completed_bounded_investigation_with_explicit_qualification_limits',
        input_hashes={a.rel(p):a.sha(p) for p in sorted(files)},
        verified_manifests=len(manifests),target_pairs=84,vela_dc_attempts=88,native_target_states=126,
        independent_column_runs=19,poisson_read_only_probes=8,
        checks={'python_regression_command':'D:/msys64/ucrt64/bin/python.exe tests/regression/test_simplemos_jacobian_fullfield_followup.py',
            'python_regression_tests':10,'python_regression_passed':10,'regression_result':'Observed in tool execution before sealing; no solver-core regression suite claimed.',
            'isolated_cpp_builds':2,'python_syntax_check':'Passed for investigation scripts','report_links':'Checked','plot_visual_check':True},
        limitations=['33 active local column/block FD signals remain unstable.',
            'Only the selected full neighborhoods were audited, not every global column or nonsmooth branch point.',
            'All shallow-channel hole rows are outside the current carrier gate; source-scale global qualification remains false.',
            'Cross-code physical current discrepancy remains; no physical correction candidate or production Jacobian promoted.'],
        production_solver_changed=False,m82_released=False,m83_released=False))
    a.verify(b.OUT/'followup_evidence.json')
    print('Sealed',len(files),'files;',len(manifests),'input manifests verified; 84 qualified pairs',flush=True)


if __name__=='__main__':main()
