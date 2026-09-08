"""Independently audit frozen gates, physical row metrics and paired comparisons, then seal."""
from pathlib import Path
import ast
import math
import re
import run_simplemos_reference_subset_20260906 as r

a=r.a;d=r.d;OUT=r.OUT;LOCAL=r.LOCAL


def main():
    manifests=[r.p.prior.OUT/'validation_evidence.json',OUT/'native_freeze.json',OUT/'export_contract.json',OUT/'vela_freeze.json',OUT/'mobility_nodal_diagnostic.json']
    for p in manifests:a.verify(p)
    jobs=a.read(OUT/'vela_contract.json')['jobs'];rows=a.rows(OUT/'vela_runs.csv');assert len(jobs)==len(rows)==24
    files=set(manifests);checks=[];qualified_counts={arm:0 for arm in r.ARMS}
    for job,reported in zip(jobs,rows):
        dest=Path(job['config']).parent;cfg=a.read(dest/'config.json');s=a.read(dest/'config.status.json');q=a.read(dest/'all_row.status.json')
        original=a.read(r.p.prior.v.m.LOCAL/job['case']/'vg_020/native/config.json')
        # Entire unchanged numerical/physical configuration, not a hand-picked subset of tolerances.
        x=cfg['solver'].copy();x.pop('local_update_diagnostics')
        y=original['solver'].copy()
        if job['model']=='masetti':y['mobility']=dict(model='masetti',doping_concentration_basis='total_impurity')
        assert x==y
        for key in ('mesh_file','node_doping_file','materials_file','scaling','doping'):assert cfg[key]==original[key]
        expect=[dict(c, bias=job['vg']) if c['name']=='gate' else c for c in original['contacts']]
        assert cfg['contacts']==expect
        e=r.env(job['arm'],dest)
        assert e['VELA_LINEAR_SOLVER']=='sparselu'
        if job['model']=='masetti':assert all(k not in e for k in ('VELA_VALIDATE_VECTOR_CHAIN_ENABLE','VELA_VALIDATE_VECTOR_CHAIN_STEP','VELA_VALIDATE_TRANSPORT_STEP'))
        assert ('VELA_SIMPLEMOS_LINEAR_REFINEMENT_PREFIX' in e)==job['arm'].endswith('_refined')
        if job['arm'].endswith('_refined'):
            linear=a.rows(dest/'linear_summary.csv');iters=sorted({int(t['iteration']) for t in linear})
            for it in iters:assert [int(t['correction']) for t in linear if int(t['iteration'])==it]==list(range(5))
        else:assert not (dest/'linear_summary.csv').exists()
        geo=d.matrix.spatial.m73.Geometry(job['device']);mask=d.matrix.spatial.old.m78.supports(job['device'],geo,.05)[0]['all_si'].copy();mask[geo.contact_nodes]=False
        terms=d.ordered(dest/'all_row.csv',geo.count);ratios=[];bad={'electron':0,'hole':0}
        for t,keep in zip(terms,mask):
            if not keep:continue
            for c in bad:
                scale=max(float(t[c+'_flux_abs_sum']),abs(float(t[c+'_recombination'])),abs(float(t[c+'_impact'])))
                assert scale>0
                ratio=abs(float(t[c+'_residual']))/scale;ratios.append(ratio);bad[c]+=ratio>1e-6
        gate=q['carrier_row_convergence'];assert len(ratios)==gate['qualified_row_count']==1814
        assert sum(bad.values())==gate['violation_count']==int(reported['all_row_violations'])
        assert all(bad[c]==int(reported[c+'_violations']) for c in bad)
        assert math.isclose(max(ratios),gate['max_ratio'],rel_tol=1e-11,abs_tol=1e-18)
        cc=s['contact_currents_A_per_um'];kcl=abs(math.fsum(cc.values()))/abs(cc['drain'])
        qualified=s['converged'] and s['exit_code']==q['exit_code']==0 and not sum(bad.values()) and q['global_continuity_closure']['satisfied'] and kcl<=1e-8
        assert qualified==(reported['comparison_qualified']=='True')
        qualified_counts[job['arm']]+=qualified
        assert math.isclose(cc['drain']/job['native_Id_A_per_um']-1,float(reported['signed_Id_error_relative']),rel_tol=1e-13)
        checks.append(dict(arm=job['arm'],case=job['case'],index=job['index'],active_rows=len(ratios),recomputed_violations=sum(bad.values()),qualified=qualified,unchanged_contract=True))
    for pair in a.rows(OUT/'nwell_pairs.csv'):
        lo=next(x for x in rows if x['arm']==pair['arm'] and x['device']=='n19' and x['vd']==pair['vd'] and x['vg']==pair['vg'])
        hi=next(x for x in rows if x['arm']==pair['arm'] and x['device']=='n23' and x['vd']==pair['vd'] and x['vg']==pair['vg'])
        # Independent expression: difference of absolute cross-TCAD log errors.
        value=math.log10(float(hi['Id_A_per_um'])/float(hi['native_Id_A_per_um']))-math.log10(float(lo['Id_A_per_um'])/float(lo['native_Id_A_per_um']))
        assert math.isclose(value,float(pair['pair_error_dex']),rel_tol=1e-12,abs_tol=1e-14)
        assert (pair['qualified']=='True')==(lo['comparison_qualified']==hi['comparison_qualified']=='True')
    assert qualified_counts=={k:v['qualified'] for k,v in a.read(OUT/'summary.json')['arms'].items()}
    assert 'All tests passed (10 assertions in 3 test cases)' in (LOCAL/'test_masetti.log').read_text()
    a.write_csv(OUT/'independent_qualification_check.csv',checks)
    for root in (LOCAL,OUT):files.update(p for p in root.rglob('*') if p.is_file() and p.name!='validation_evidence.json')
    scripts=[d.REPO/'scripts'/n for n in ('prepare_simplemos_reference_subset_20260906.py','run_simplemos_reference_subset_20260906.py','analyze_simplemos_reference_subset_20260906.py','audit_simplemos_subset_mobility_export_20260907.py','seal_simplemos_reference_subset_20260907.py')]
    for p in scripts:ast.parse(p.read_text(encoding='utf-8-sig'),filename=str(p))
    files.update(scripts);files.add(r.RUNNER);files.add(d.REPO/'build-release/test_mobility.exe');files.add(d.REPO/'build-release/sentaurus_import.exe')
    report=d.REPO/'docs/validation/simplemos_reference_subset_comparison_2026-09-06.md';files.add(report)
    for target in re.findall(r'\]\(([^)]+)\)',report.read_text(encoding='utf-8')):
        if '://' in target:continue
        path=(report.parent/target.split('#')[0]).resolve()
        if path.name!='validation_evidence.json':assert path.exists(),path;files.add(path)
    a.write(OUT/'validation_evidence.json',dict(schema_version=1,completed_date='2026-09-07',status='completed_bounded_reference_subset_comparison',
        input_hashes={a.rel(p):a.sha(p) for p in sorted(files)},verified_manifests=len(manifests),
        native_states=16,vela_states=24,post_row_probes=24,independent_physical_rows_checked=24*1814,qualified_counts=qualified_counts,
        checks=dict(frozen_inputs='passed',unchanged_gates='passed',physical_row_recalculation='passed',pair_identity='passed',
            native_hdf5_exports=16,masetti_catch2_tests=3,masetti_catch2_assertions=10,python_syntax=len(scripts),report_links='passed',full_ctest='not_run_no_core_change'),
        production_changes=False,acceptance_changes=False,m82_released=False,m83_released=False,
        limitations=['Eight independent working points; not a full curve or initialization-invariance qualification.',
            'Composite mobility simplification does not identify each removed model contribution.',
            'All eight no-refinement Masetti outputs remain failed diagnostics.',
            'Native nodal mobility and local pointwise formula are not proven equivalent internal quantities.']))
    a.verify(OUT/'validation_evidence.json');print('Sealed',len(files),'files; independent physical rows checked:',24*1814,flush=True)


if __name__=='__main__':main()
