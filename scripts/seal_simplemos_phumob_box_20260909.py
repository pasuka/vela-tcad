"""Summarize paired controls and preserve the PhuMob box candidate evidence."""
import argparse
import ast
import re
import shutil
import xml.etree.ElementTree as ET
from pathlib import Path
import validate_simplemos_phumob_box_qualified_20260909 as stage

a,d=stage.a,stage.d
OUT,LOCAL,REPO=stage.OUT,stage.LOCAL,stage.v.REPO
REPORT=REPO/'docs/validation/simplemos_phumob_box_candidate_validation_2026-09-09.md'
SCRIPTS=[REPO/'scripts'/name for name in ('validate_simplemos_phumob_box_20260909.py',
    'analyze_simplemos_phumob_box_20260909.py','validate_simplemos_phumob_box_qualified_20260909.py',
    'seal_simplemos_phumob_box_20260909.py')]


def review():
    a.verify(OUT/'dc_comparison_evidence.json');a.verify(OUT/'post_evidence.json')
    xml=LOCAL/'ctest_full.xml';tree=ET.parse(xml).getroot();tests=list(tree.iter('testcase'))
    failed=[t.attrib['name'] for t in tests if t.find('failure') is not None]
    old=stage.v.previous.OUT/'ctest_comparison.json';previous=a.read(old)
    log=(LOCAL/'ctest_full.log').read_text(encoding='utf8')
    seconds=float(re.search(r'Total Test time \(real\) =\s*([\d.]+)',log)[1])
    ctest=dict(total=len(tests),passed=len(tests)-len(failed),failed=failed,
        new_failures=sorted(set(failed)-set(previous['failed'])),
        missing_previous_failures=sorted(set(previous['failed'])-set(failed)),seconds=seconds,
        phumob_box_tests=sum(t.attrib['name'].startswith('PhuMob box ') for t in tests))
    assert not ctest['new_failures'] and not ctest['missing_previous_failures'],ctest
    assert ctest['phumob_box_tests']==3
    a.write(OUT/'ctest_comparison.json',ctest)
    rows=a.rows(OUT/'dc_comparison.csv');pairs=[]
    for vd,vg in sorted({(r['vd'],r['vg']) for r in rows}):
        lo=next(r for r in rows if r['vd']==vd and r['vg']==vg and r['device']=='n19')
        hi=next(r for r in rows if r['vd']==vd and r['vg']==vg and r['device']=='n23')
        native_ratio=float(hi['native_Id_A_per_um'])/float(lo['native_Id_A_per_um'])
        native_delta=float(hi['native_Id_A_per_um'])-float(lo['native_Id_A_per_um'])
        for arm in ('legacy','candidate'):
            ratio=float(hi[arm+'_Id_A_per_um'])/float(lo[arm+'_Id_A_per_um'])
            delta=float(hi[arm+'_Id_A_per_um'])-float(lo[arm+'_Id_A_per_um'])
            pairs.append(dict(vd=vd,vg=vg,arm=arm,native_high_over_low=native_ratio,
                high_over_low=ratio,ratio_error_percent=100*(ratio/native_ratio-1),
                native_high_minus_low_A_per_um=native_delta,high_minus_low_A_per_um=delta,
                delta_error_A_per_um=delta-native_delta))
    a.write_csv(OUT/'nwell_pairs.csv',pairs)
    attempts=[r for arm in ('legacy','candidate') for r in a.rows(OUT/arm/'attempts.csv')]
    accepted=[r for r in attempts if r['qualified']=='True']
    preflight=[a.read(p) for p in (LOCAL/'fixed').glob('*/density_preflight.json')]
    facts=dict(points=len(rows),attempts=len(attempts),accepted_attempts=len(accepted),
        failed_attempts=len(attempts)-len(accepted),max_row_ratio=max(float(r['max_row_ratio']) for r in accepted),
        max_kcl_over_Id=max(float(r['kcl_over_Id']) for r in accepted),
        max_port_relative=max(float(r['port_relative']) for r in accepted),
        fixed_state_density_preflight_count=len(preflight),
        fixed_state_density_max_relative=max(r['max_density_relative'] for r in preflight),
        distinct_nonzero_cross_entries=len({(r['key'],r['direction'].rsplit('_',1)[0],r['output_block'],r['output_node']) for r in a.rows(OUT/'cross_checks.csv')}),
        native_error_improved_points=sum(abs(float(r['candidate_versus_native_percent']))<abs(float(r['legacy_versus_native_percent'])) for r in rows),
        candidate_pair_error_max_percent=max(abs(r['ratio_error_percent']) for r in pairs if r['arm']=='candidate'),
        legacy_pair_error_max_percent=max(abs(r['ratio_error_percent']) for r in pairs if r['arm']=='legacy'))
    a.write(OUT/'review_facts.json',facts)
    d.matrix.freeze(OUT/'review_evidence.json',[OUT/'dc_comparison_evidence.json',OUT/'post_evidence.json',
        xml,LOCAL/'ctest_full.log',old,OUT/'ctest_comparison.json',OUT/'nwell_pairs.csv',OUT/'review_facts.json',Path(__file__).resolve()])
    print(ctest,facts,pairs,flush=True)


def seal():
    for script in SCRIPTS:ast.parse(script.read_text(encoding='utf8'))
    identities={};manifests=[]
    for file in stage.PARENT_OUT.rglob('*.json'):
        data=a.read(file)
        if isinstance(data,dict) and 'input_hashes' in data:
            manifests.append(file)
            for name,expected in data['input_hashes'].items():
                assert name not in identities or identities[name]==expected,(name,file)
                identities[name]=expected
    for name,expected in identities.items():assert a.sha(REPO/name)==expected,name
    links=re.findall(r'\]\(([^)]+)\)',REPORT.read_text(encoding='utf8'))
    for link in links:assert (REPORT.parent/link).resolve().exists(),link
    production=[REPO/name for name in ('include/vela/physics/MobilityModel.h','src/physics/MobilityModel.cpp',
        'include/vela/equation/AssemblerUtils.h','include/vela/equation/CoupledDDAssembler.h',
        'src/equation/CoupledDDAssembler.cpp','src/equation/FixedStateOperatorAudit.cpp',
        'src/post/ContactCurrent.cpp','src/tools/vela_example_runner.cpp','tests/test_element_box_transport.cpp',
        'docs/config_schema.md','docs/validation/simplemos_branch_status.md')]
    extra=[REPORT,stage.q.run.RUNNER,REPO/'build-release/libvela_core.a',REPO/'build-release/CMakeCache.txt',
        OUT/'review_evidence.json',OUT/'dc_summary.json',OUT/'fixed_summary.json',OUT/'post_summary.json',
        OUT/'initial_state_failure_ledger.csv',stage.PARENT_OUT/'fixed_summary.json',stage.PARENT_OUT/'probe_summary.json']
    root=REPO/'build-release/phumob_box_snapshot_20260909';root.mkdir(exist_ok=False);snapshots=[]
    for i,file in enumerate(SCRIPTS+production+extra):
        target=root/f'{i:02d}_{file.name}';shutil.copy2(file,target);snapshots.append(target)
    quality=dict(python_syntax_checks=len(SCRIPTS),manifest_identity_checks=len(manifests),
        unique_file_hashes_verified=len(identities),report_links_checked=len(links),snapshots=len(snapshots),
        explicit_candidate='element_box_phumob',global_default_changed=False,native_G_floor_matched=False,
        initial_failed_state_audit_preserved=True,full_ctest_all_green=False,commit_or_push=False)
    a.write(OUT/'completion_review.json',quality)
    d.matrix.freeze(OUT/'completion_evidence.json',manifests+SCRIPTS+production+extra+snapshots+[OUT/'completion_review.json'])
    print(quality,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('review','seal'))
    globals()[parser.parse_args().action]()
