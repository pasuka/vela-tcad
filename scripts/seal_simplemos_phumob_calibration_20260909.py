"""Recompute scale context and seal the PhuMob calibration report's evidence."""
import ast,math,re,shutil
from pathlib import Path
import analyze_simplemos_phumob_assembled_cross_20260909 as cross
import analyze_simplemos_phumob_native_cells_20260909 as cell

a,d=cross.a,cross.d;p=cross.p


def main():
    sources=[cell.s.OUT/'native_evidence.json',cell.s.OUT/'export_evidence.json',cell.OUT/'evidence.json',
             cell.OUT/'screening_floor/evidence.json',cross.OUT/'evidence.json',p.OUT/'cross_derivatives/runtime_ready/evidence.json']
    for source in sources:a.verify(source)
    rows=a.rows(cross.OUT/'cross_entries.csv')
    unique={(r['case'],r['input_mode'],r['input_node'],r['output_block'],r['output_node']):r for r in rows}
    assert len(rows)==328 and len(unique)==164
    assert sum(r['qualified']=='False' for r in unique.values())==106
    assert sum(r['production_zero']=='True' for r in unique.values())==94
    context=[]
    for job in a.read(cross.original.OUT/'contract.json')['jobs']:
        sampled=a.rows(Path(job['dir'])/'rows.csv')
        for node in cross.NODES:
            for mode in ('phin','phip'):
                direct=math.fsum(abs(float(r['analytic_derivative'])) for r in sampled if r['direction']==f'{mode}_{node}_1e-05' and r['row_block']==mode)
                subset=[r for r in unique.values() if r['case']==job['case'] and r['input_mode']==mode and int(r['input_node'])==node]
                error=math.fsum(abs(float(r['absolute_error'])) for r in subset)
                ref=math.fsum(abs(float(r['reference'])) for r in subset)
                context.append(dict(case=job['case'],input_mode=mode,input_node=node,sampled_direct_column_L1=direct,cross_error_L1=error,cross_reference_L1=ref,error_over_direct=error/direct if direct else None))
    context_path=cross.OUT/'column_scale_context.csv'
    if context_path.exists():
        previous=a.rows(context_path)
        assert len(previous)==len(context)
        for old,new in zip(previous,context):
            assert all(old[k]==('' if v is None else str(v)) for k,v in new.items())
    else:a.write_csv(context_path,context)
    report=p.REPO/'docs/validation/simplemos_phumob_cell_edge_and_cross_calibration_2026-09-09.md'
    status=p.REPO/'docs/validation/simplemos_branch_status.md'
    files=[p.REPO/'scripts'/name for name in ('simplemos_phumob_calibration_native_20260908.py','simplemos_phumob_supported_export_20260908.py','audit_simplemos_phumob_cross_derivatives_20260908.py',
        'analyze_simplemos_phumob_native_cells_20260909.py','audit_simplemos_phumob_assembled_cross_20260909.py','analyze_simplemos_phumob_assembled_cross_20260909.py',
        'localize_simplemos_phumob_screening_floor_20260909.py','seal_simplemos_phumob_calibration_20260909.py')]
    for file in files:ast.parse(file.read_text(encoding='utf8'))
    links=re.findall(r'\]\(([^)]+)\)',report.read_text(encoding='utf8'))
    for link in links:
        assert (report.parent/link).resolve().exists(),link
    assert all(r['native_qualified']=='True' for r in a.rows(cell.s.OUT/'native_points.csv'))
    summary=a.read(cell.OUT/'summary.json')
    assert summary['chosen_cell_failures']==8 and summary['port_replay_failures']==0
    floor=a.read(cell.OUT/'screening_floor/summary.json')
    assert floor['native_floor_algorithm_identified'] is False
    assert floor['heldout_diagnostic_max_relative']<1e-7
    quality=dict(python_syntax_checks=len(files),evidence_identity_checks=len(sources),report_links_checked=len(links),
        native_states_qualified=16,port_replays_qualified=16,production_element_box_phumob_qualified=False,
        matrix_distinct_entries=164,matrix_distinct_failures=106,matrix_distinct_zeros=94,
        max_cross_error_over_sampled_direct_column=max(r['error_over_direct'] for r in context if r['error_over_direct'] is not None),
        interpretation='Completion of calibration and localization; failed numerical gates remain failed. No new production solver build, full CTest, PhuMob self-consistent solve, commit, or push.')
    quality_file=p.OUT/'completion_review_20260909.json'
    if quality_file.exists():assert a.read(quality_file)==quality
    else:a.write(quality_file,quality)
    extra=[report,status,p.REPO/'scripts/diagnostics/simplemos_phumob_scalar_probe.cpp',cross.OUT/'column_scale_context.csv',quality_file,
           cell.s.LOCAL/'complete.txt',cell.s.LOCAL/'guard_process.json',p.OUT/'initial_attempt_ledger.csv',
           cell.OUT/'preflight_failure/failure.json',cross.original.OUT/'analysis_preflight_failure.json']
    snapshots=[]
    for index,file in enumerate(files+extra):
        # Keep snapshots shallow to stay within Windows path-length limits.
        dest=p.REPO/'build-release/phumob_snapshot_20260909'/f'{index:02d}_{file.name}'
        dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(file,dest);snapshots.append(dest)
    d.matrix.freeze(p.OUT/'completion_evidence_20260909.json',sources+files+extra+snapshots)
    print(quality,flush=True)


if __name__=='__main__':main()
