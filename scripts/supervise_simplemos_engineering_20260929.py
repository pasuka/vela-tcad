"""Continue the sealed engineering matrix after its first curve passes.

Continuous compute job, not a scheduler. Two workers at most; failures stop
dispatch without killing another in-flight solver or changing any gate.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
import ctypes
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import simplemos_hfs_cloud_20260926 as h


def wait_process(pid):
    if not pid:
        return
    if sys.platform != 'win32':
        raise ValueError('The supplied pilot PID must belong to this Windows host')
    api=ctypes.WinDLL('kernel32',use_last_error=True)
    api.OpenProcess.argtypes=[ctypes.c_ulong,ctypes.c_int,ctypes.c_ulong]
    api.OpenProcess.restype=ctypes.c_void_p
    api.WaitForSingleObject.argtypes=[ctypes.c_void_p,ctypes.c_ulong]
    api.WaitForSingleObject.restype=ctypes.c_ulong
    api.CloseHandle.argtypes=[ctypes.c_void_p]
    handle=api.OpenProcess(0x00100000,False,pid)
    if not handle:
        if ctypes.get_last_error()==87:  # Already exited; summary still required.
            return
        raise ctypes.WinError(ctypes.get_last_error())
    try:
        while True:
            result=api.WaitForSingleObject(handle,10000)
            if result==0:
                return
            if result!=258:
                raise ctypes.WinError(ctypes.get_last_error())
    finally:
        api.CloseHandle(handle)


def verified_summary(base,case):
    path=base/'matrix'/case/'summary.json'
    if not path.exists():
        raise ValueError(f'No completed summary for {case}')
    data=h.read(path)
    rows=data.get('results',[])
    if not data.get('passed') or data.get('points')!=51 or len(rows)!=51:
        raise ValueError(f'Incomplete or failed case {case}')
    if {r['index'] for r in rows}!=set(range(51)) or any(r['case']!=case or not r['passed'] for r in rows):
        raise ValueError(f'Invalid case coverage {case}')
    return rows


def run(base,pilot_pid):
    base=base.resolve();root=base/'supervisor';root.mkdir(exist_ok=True)
    native=h.read(base/'native_verified.json')
    if not native.get('passed') or native.get('runs')!=16 or native.get('states')!=816:
        raise ValueError('Complete native archive verification is required')
    # Never silently reuse a stale lock after an interrupted run.
    with (root/'execution.lock').open('x') as f:
        json.dump(dict(pid=os.getpid(),started=time.time()),f)
    driver=base/'source/scripts/run_simplemos_engineering_matrix_20260929.py'
    source=base/'source';runner=source/'build-release/vela_example_runner.exe'
    cases=[c['case'] for c in h.read(base/'data/inputs/contract.json')['cases']]
    if len(cases)!=16 or len(set(cases))!=16:
        raise ValueError('Expected sixteen unique cases')
    watched=[driver,Path(__file__),runner,*[source/'scripts'/name for name in (
        'simplemos_hfs_cloud_20260926.py','simplemos_original_matrix_cloud_20260926.py',
        'check_simplemos_cold_overlap_20260927.py','simplemos_srh_cloud_20260926.py',
        'sentaurus_import.py','state_archive.py')]]
    hashes={str(p):h.sha(p) for p in watched}
    h.write(root/'seal.json',dict(hashes=hashes,cases=cases,pilot_pid=pilot_pid,max_workers=2))
    status=dict(phase='waiting_for_pilot',completed=[],running=[],failed=[],requested_points=816)
    h.write(root/'status.json',status)
    wait_process(pilot_pid)
    pilot='n17_vd_0p05'
    results=verified_summary(base,pilot);status['completed'].append(pilot)

    def one(case):
        for path,expected in hashes.items():
            if h.sha(path)!=expected:
                raise ValueError(f'Sealed executable or helper changed: {path}')
        with (root/(case+'.stdout.log')).open('x') as out,(root/(case+'.stderr.log')).open('x') as err:
            result=subprocess.run([sys.executable,'-B','-X','utf8',str(driver),'--base',str(base),'--case',case],cwd=source,stdout=out,stderr=err)
        if result.returncode:
            raise ValueError(f'{case} exited {result.returncode}; inspect preserved case evidence')
        return verified_summary(base,case)

    pending=iter(c for c in cases if c!=pilot);futures={};status['phase']='running'
    with ThreadPoolExecutor(max_workers=2) as pool:
        def dispatch():
            while len(futures)<2 and not status['failed']:
                case=next(pending,None)
                if case is None:
                    return
                futures[pool.submit(one,case)]=case
        dispatch()
        while futures:
            status['running']=list(futures.values());h.write(root/'status.json',status)
            finished,_=wait(futures,timeout=30,return_when=FIRST_COMPLETED)
            for future in finished:
                case=futures.pop(future)
                try:
                    results.extend(future.result());status['completed'].append(case)
                except Exception as exc:
                    status['failed'].append(dict(case=case,error=repr(exc)))
            dispatch()
    status.update(phase='failed' if status['failed'] else 'completed',running=[])
    h.write(root/'status.json',status)
    if status['failed']:
        raise ValueError('Matrix stopped on a preserved failure; no gates were changed')
    if len(results)!=816 or len({(r['case'],r['index']) for r in results})!=816:
        raise ValueError('Final matrix has incomplete or duplicate coverage')
    field_rows=[]
    for case in cases:
        for i in range(51):
            point=h.read(base/'matrix'/case/f'vg_{i:03d}/fields.json')
            field_rows.extend(dict(case=case,index=i,**r) for r in point['fields'])
            field_rows.append(dict(case=case,index=i,field='SRH',unit='weighted_relative_L1',max_abs=point['SRH_weighted_L1']))
    h.csvout(root/'physical_fields.csv',field_rows)
    h.csvout(root/'comparison.csv',[dict(case=r['case'],index=r['index'],vg=r['vg'],
        cold_Id_error_percent=r['Id_error_percent']['cold'],native_Id_error_percent=r['Id_error_percent']['native'],
        supplemental_native_Id_change_percent=r['supplemental_native_Id_change_percent'],**r['dual']) for r in results])
    h.write(root/'summary.json',dict(points=816,passed=True,
        max_abs_Id_error_percent=max(abs(v) for r in results for v in r['Id_error_percent'].values()),
        max_dual_Id_relative=max(r['dual']['Id_relative'] for r in results),
        remaining='Full regression and review/commit precede performance work; mobility field definition coverage is separate'))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--base',type=Path,required=True);p.add_argument('--pilot-pid',type=int,default=0);a=p.parse_args()
    try:
        run(a.base,a.pilot_pid)
    except Exception as exc:
        h.write(a.base/'supervisor/failure.json',dict(error=repr(exc),time=time.time()))
        raise
