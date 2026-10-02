"""Gated split continuation requalification on the existing Codespace snapshot."""
import argparse,json,re,time,subprocess,sys
from pathlib import Path
import simplemos_hfs_cloud_20260926 as h
from simplemos_srh_cloud_20260926 import Experiment
from analyze_simplemos_srh_cloud_20260926 import analyze
from simplemos_original_matrix_cloud_20260926 import Matrix
from check_simplemos_cold_overlap_20260927 import compare as compare_cold
from check_simplemos_continuation_identity_20260927 import compare as compare_identity

def run(base,revision=''):
    source=base/('source-continuation'+revision);build=source/'build'
    root=base/('continuation_validation_20260927'+revision);root.mkdir(exist_ok=True)
    build_name='build_continuation'+revision
    test_name='continuation'+revision+'_ctest'
    def progress(stage,**more):
        h.write(root/'progress.json',dict(stage=stage,**more));print(stage,more,flush=True)
    progress('waiting_build_and_regression')
    while not (base/(build_name+'.exit')).exists():time.sleep(10)
    assert (base/(build_name+'.exit')).read_text().strip()=='0','Build or split runtime tests failed'
    while not (base/(test_name+'.exit')).exists():time.sleep(10)
    failures=base/('continuation'+revision+'_full_failures.txt')
    if not failures.exists():failures=build/'Testing/Temporary/LastTestsFailed.log'
    actual=sorted(line.split(':',1)[1] for line in failures.read_text().splitlines()) if failures.exists() else []
    resolved=[]
    retest=base/('continuation'+revision+'_environment_retest.log')
    if retest.exists():
        report=retest.read_text()
        assert '1/1 Test' in report and 'pn2d_bv_process_observability_regression' in report
        assert '100% tests passed, 0 tests failed out of 1' in report
        for row in h.read(base/'continuation_restored_files.json'):
            assert h.sha(source/row['path'])==row['sha256']
        resolved=['pn2d_bv_process_observability_regression']
    old=h.read(base/'srh_regression.json')['failures'];new=sorted(set(actual)-set(old)-set(resolved))
    log=(base/(test_name+'.log')).read_text()
    total=int(re.search(r'out of (\d+)',log).group(1))
    regression=dict(total=total,initial_failures=actual,failures=sorted(set(actual)-set(resolved)),resolved_environment_failures=resolved,new_failures=new)
    h.write(root/'regression.json',regression)
    assert total>=984 and not new,'Regression has new failures'
    runner=build/'vela_example_runner'
    progress('sixteen_point_dual',regression=regression)
    experiment=Experiment(root/'controls',base/'replay_v4',runner)
    experiment.finite();analyze(root/'controls',base/'replay_v4')
    compare_identity(base,revision)
    progress('cold_start_pilot',control_summary=h.read(root/'controls/summary.json'))
    pilot=Matrix(root/'original_pilot',base/'original_inputs',runner,root/'controls')
    pilot.all(True)
    compare_cold(root/'controls',root/'original_pilot',base/'original_inputs')
    progress('original_816',pilot_summary=h.read(root/'original_pilot/summary.json'))
    matrix=Matrix(root/'original_matrix',base/'original_inputs',runner,root/'controls')
    matrix.all(False)
    compare_cold(root/'controls',root/'original_matrix',base/'original_inputs')
    progress('complete',matrix_summary=h.read(root/'original_matrix/summary.json'))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--base',type=Path,required=True);p.add_argument('--revision',default='');a=p.parse_args()
    try:run(a.base.resolve(),a.revision)
    except Exception as e:
        root=a.base/('continuation_validation_20260927'+a.revision);root.mkdir(exist_ok=True)
        h.write(root/'failure.json',dict(error=repr(e),progress=h.read(root/'progress.json') if (root/'progress.json').exists() else None))
        raise
