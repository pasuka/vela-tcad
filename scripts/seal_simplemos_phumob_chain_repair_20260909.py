"""Review and preserve the scoped PhuMob chain repair, including failed controls."""
import argparse
import ast
import re
import shutil
import xml.etree.ElementTree as ET
from pathlib import Path

import validate_simplemos_phumob_chain_fix_20260909 as v
import simplemos_phumob_floor_ramp_20260909 as ramp

a, d, OUT, LOCAL, REPO = v.a, v.d, v.OUT, v.LOCAL, v.p.REPO
REPORT = REPO / 'docs/validation/simplemos_phumob_stable_chain_repair_2026-09-09.md'
STATUS = REPO / 'docs/validation/simplemos_branch_status.md'
SCRIPTS = [REPO / 'scripts' / name for name in (
    'validate_simplemos_phumob_chain_fix_20260909.py',
    'analyze_simplemos_phumob_chain_fix_20260909.py',
    'simplemos_phumob_floor_temperature_20260909.py',
    'simplemos_phumob_floor_ramp_20260909.py',
    'analyze_simplemos_phumob_floor_temperature_20260909.py',
    'audit_simplemos_phumob_native_cutoff_20260909.py',
    'seal_simplemos_phumob_chain_repair_20260909.py')]


def review():
    # CTest resolves --output-junit relative to its preset build directory.
    xml = REPO / 'build-release/build-release/phumob_fix_20260909/ctest_full.xml'
    tree = ET.parse(xml).getroot()
    tests = list(tree.iter('testcase'))
    failed = [t for t in tests if t.find('failure') is not None]
    old_file = REPO / 'reference_tcad/simplemos_sentaurus2022/bgn_psi_derivative_20260908/state_independent/ctest_comparison.json'
    old = a.read(old_file)
    names = [t.attrib['name'] for t in failed]
    categories = []
    for test in failed:
        output = test.findtext('system-out', '')
        if 'FileNotFoundError' in output and 'test_mos_mixed_material.cpp' in output:
            category = 'historical_missing_fixture'
        elif 'not frozen by' in output or ('assertEqual(expected, sha256' in output):
            category = 'historical_identity_or_supersession'
        else:
            raise AssertionError((test.attrib['name'], output))
        categories.append(dict(test=test.attrib['name'], category=category))
    comparison = dict(total=len(tests), passed=len(tests)-len(failed), failed=names,
        new_failures=sorted(set(names)-set(old['failed'])),
        missing_previous_failures=sorted(set(old['failed'])-set(names)),
        phumob_tests=sum('PhuMob' in t.attrib['name'] for t in tests),
        phumob_failed=sum('PhuMob' in t.attrib['name'] for t in failed),
        seconds=253.10, junit_path=a.rel(xml))
    assert comparison['total'] == 785 and comparison['passed'] == 769
    assert not comparison['new_failures'] and not comparison['missing_previous_failures']
    assert comparison['phumob_tests'] == 18 and comparison['phumob_failed'] == 0
    a.write(OUT/'ctest_comparison.json', comparison)
    a.write_csv(OUT/'ctest_failure_categories.csv', categories)

    entries = a.rows(OUT/'cross_analysis/cross_entries.csv')
    distinct = {(r['case'], r['input_mode'], r['input_node'], r['output_block'], r['output_node']): r for r in entries}
    assert len(entries) == 328 and len(distinct) == 164
    assert all(r['qualified'] == 'True' and r['production_zero'] == 'False' for r in entries)
    before, after = [a.rows(OUT/arm/'attempts.csv') for arm in ('before', 'after')]
    assert len(before) == len(after) == 8
    assert all(r['qualified'] == 'True' and r['attempt'] == '0' for r in before+after)
    assert [r['iterations'] for r in before] == [r['iterations'] for r in after]
    cutoff = a.rows(ramp.OUT/'cutoff/results.csv')
    for key in ('T299', 'T300'):
        selected = [r for r in cutoff if r['key'] == key and r['carrier'] == 'h']
        assert [r['mode'] for r in selected if r['distinction_gate_passed'] == 'True'] == ['P_lower']
    selected = [r for r in cutoff if r['key'] == 'T350' and r['carrier'] == 'h']
    assert {r['mode'] for r in selected if r['distinction_gate_passed'] == 'True'} == {'P_upper', 'G_floor'}
    facts = dict(distinct_cross_entries=len(distinct), distinct_cross_failures=0,
        cross_max_relative=max(float(r['relative_error']) for r in entries),
        before_after_max_row_ratio=max(float(r['max_row_ratio']) for r in before+after),
        before_after_max_kcl_over_Id=max(float(r['kcl_over_Id']) for r in before+after),
        before_after_max_port_relative=max(float(r['port_relative']) for r in before+after),
        iterations_unchanged=True, native_internal_minimizer_identified=False,
        native_cutoff_result='At 299/300 K holes select the lower root; at 350 K holes exclude the lower root but cannot distinguish the upper root from G-value clipping. Electron samples do not distinguish the hypotheses.')
    a.write(OUT/'review_facts.json', facts)
    d.matrix.freeze(OUT/'review_evidence.json', [xml, LOCAL/'ctest_full.log', old_file,
        OUT/'ctest_comparison.json', OUT/'ctest_failure_categories.csv', OUT/'review_facts.json',
        OUT/'cross_analysis/cross_entries.csv', OUT/'before/attempts.csv', OUT/'after/attempts.csv',
        ramp.OUT/'cutoff/results.csv', Path(__file__).resolve()])
    print(comparison, facts, flush=True)


def seal():
    for script in SCRIPTS:
        ast.parse(script.read_text(encoding='utf8'))
    manifests = []
    identities = {}
    for root in (OUT, ramp.OUT.parent):
        for file in root.rglob('*.json'):
            data = a.read(file)
            if isinstance(data, dict) and 'input_hashes' in data:
                manifests.append(file)
                for name, expected in data['input_hashes'].items():
                    assert name not in identities or identities[name] == expected, (name, file)
                    identities[name] = expected
    for name, expected in identities.items():
        assert a.sha(REPO/name) == expected, name
    links = re.findall(r'\]\(([^)]+)\)', REPORT.read_text(encoding='utf8'))
    for link in links:
        assert (REPORT.parent/link).resolve().exists(), link
    production = [REPO/name for name in (
        'include/vela/physics/MobilityModel.h', 'src/physics/MobilityModel.cpp',
        'include/vela/equation/AssemblerUtils.h', 'src/equation/CoupledDDAssembler.cpp',
        'include/vela/discretization/StableSGDerivative.h', 'tests/test_phumob.cpp',
        'docs/config_schema.md')]
    extras = [REPORT, STATUS, OUT/'review_evidence.json', OUT/'repair_summary.json',
        ramp.OUT.parent/'temperature_preflight_failure.csv',
        ramp.OUT.parent/'upload_scope_review.json', ramp.OUT/'analysis/source_identity.json',
        v.RUNNER, REPO/'build-release/libvela_core.a', REPO/'build-release/CMakeCache.txt']
    snapshot_root = REPO/'build-release/phumob_fix_snapshot_20260909'
    snapshot_root.mkdir(exist_ok=False)
    snapshots = []
    for i, file in enumerate(SCRIPTS+production+extras):
        target = snapshot_root/f'{i:02d}_{file.name}'
        shutil.copy2(file, target)
        snapshots.append(target)
    quality = dict(python_syntax_checks=len(SCRIPTS), evidence_manifests=len(manifests),
        unique_file_hashes_verified=len(identities), report_links_checked=len(links),
        snapshot_files=len(snapshots), production_jacobian_repaired=True,
        native_cutoff_internal_algorithm_identified=False,
        native_element_box_phumob_restored=False, full_ctest_all_green=False,
        commit_or_push=False)
    a.write(OUT/'completion_review.json', quality)
    d.matrix.freeze(OUT/'completion_evidence.json', manifests+SCRIPTS+production+extras+snapshots+[OUT/'completion_review.json'])
    print(quality, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=('review', 'seal'))
    globals()[parser.parse_args().action]()
