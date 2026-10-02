"""Sealed Linux dispatch of the thirteen cases never started on T470p.

This queue stops dispatching after a failure, drains in-flight work, and never
interprets its partial assignment as completion of the sixteen-case matrix.
"""
import argparse
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import time

import simplemos_hfs_cloud_20260926 as h
from supervise_simplemos_engineering_20260929 import verified_summary


def validate_assignment(a,contract):
    all_cases={c['case'] for c in contract['cases']}
    cases=a['cases'];retained=set(a['retained_T470p'])
    if len(cases)!=13 or len(set(cases))!=13 or retained!={'n17_vd_0p05','n17_vd_1','n18_vd_0p05'}:
        raise ValueError('Unexpected ownership or duplicate dispatch')
    if set(cases)&retained or set(cases)|retained!=all_cases or a['points']!=663 or a['max_workers']!=2:
        raise ValueError('Assignment does not partition the frozen matrix')
    return cases


def checked_path(root,name):
    path=(root/name).resolve()
    if not path.is_relative_to(root.resolve()):raise ValueError('Path outside evidence root')
    return path


def verify_hashes(base,manifest):
    for rel,expected in manifest.items():
        if h.sha(checked_path(base,rel))!=expected:raise ValueError('Hash mismatch: '+rel)


def verify_source(base,a):
    old=Path(a['frozen_source']);manifest=h.read(base/'core_source_identity.json')
    actual={p.relative_to(old).as_posix() for name in ('src','include') for p in (old/name).rglob('*') if p.is_file()}
    if actual!=set(manifest):raise ValueError('Frozen core source file coverage differs')
    raw=0
    for rel,expected in manifest.items():
        data=(old/rel).read_bytes()
        raw+=hashlib.sha256(data).hexdigest()==expected['raw']
        if hashlib.sha256(data.replace(b'\r\n',b'\n')).hexdigest()!=expected['lf']:
            raise ValueError('Substantive core source difference: '+rel)
    runner=old/'build/vela_example_runner';importer=old/'build/sentaurus_import'
    if h.sha(runner)!=a['runner_sha256'] or h.sha(importer)!=a['importer_sha256']:
        raise ValueError('Frozen binary identity mismatch')
    h.write(base/'source_verified.json',dict(passed=True,files=len(manifest),raw_identical=raw,
        normalized_identical=len(manifest),runner=str(runner),runner_sha256=h.sha(runner),
        importer=str(importer),importer_sha256=h.sha(importer),
        build='Existing frozen GCC16 Release build; no solver rebuild or numerical change'))
    return runner,importer


def verify_native(base,a):
    root=base/'native';root.mkdir(exist_ok=True);files=0
    for name,expected in h.read(base/'native_archives.json').items():
        archive=base/'archives'/name
        if h.sha(archive)!=expected:raise ValueError('Native archive mismatch: '+name)
        with tarfile.open(archive) as tar:
            for member in tar.getmembers():
                checked_path(root,member.name)
                if not member.isfile() and not member.isdir():raise ValueError('Unsupported archive member')
            tar.extractall(root,filter='data')
        for line in (root/(Path(name).stem+'.sha256')).read_text().splitlines():
            expected,rel=line.split(maxsplit=1);rel=rel.lstrip('*')
            if h.sha(checked_path(root,rel))!=expected:raise ValueError('Native file mismatch: '+rel)
            files+=1
    for case in a['cases']:
        prefix='m60fields_'+case
        if (root/(prefix+'.exitcode')).read_text().strip()!='0':raise ValueError('Native exit not zero: '+case)
        states=list((root/'bundle'/case.split('_')[0]).glob(prefix+'_state_*_des.tdr'))
        if len(states)!=51:raise ValueError('Missing native states: '+case)
    h.write(base/'native_verified.json',dict(passed=True,files=files,assigned_curves=13,assigned_states=663,
        archive_curves=14,extra_curve='n18_vd_0p05 input only; not dispatched'))


def collect(base,cases):
    results=[];fields=[]
    for case in cases:
        results.extend(verified_summary(base,case))
        for index in range(51):
            point=h.read(base/'matrix'/case/f'vg_{index:03d}/fields.json')
            fields.extend(dict(case=case,index=index,**r) for r in point['fields'])
            fields.append(dict(case=case,index=index,field='SRH',unit='weighted_relative_L1',max_abs=point['SRH_weighted_L1']))
    if len(results)!=663 or len({(r['case'],r['index']) for r in results})!=663:
        raise ValueError('Incomplete cloud coverage')
    h.csvout(base/'supervisor/physical_fields.csv',fields)
    h.csvout(base/'supervisor/comparison.csv',[dict(case=r['case'],index=r['index'],vg=r['vg'],
        cold_Id_error_percent=r['Id_error_percent']['cold'],native_Id_error_percent=r['Id_error_percent']['native'],
        supplemental_native_Id_change_percent=r['supplemental_native_Id_change_percent'],**r['dual']) for r in results])
    h.write(base/'supervisor/summary.json',dict(passed=True,points=663,curves=13,
        full_816_completed=False,max_abs_Id_error_percent=max(abs(v) for r in results for v in r['Id_error_percent'].values()),
        max_dual_Id_relative=max(r['dual']['Id_relative'] for r in results),
        remaining='Join T470p results and resolve its n17 high-Vd runtime failure before full regression, review/commit and performance'))


def run(base):
    base=base.resolve();root=base/'supervisor';root.mkdir(exist_ok=True)
    with (root/'execution.lock').open('x') as f:f.write(str(os.getpid()))
    a=h.read(base/'assignment.json');cases=validate_assignment(a,h.read(base/'data/inputs/contract.json'))
    state=dict(phase='verifying',pid=os.getpid(),queued=cases.copy(),running={},completed=[],failed=[],
               requested_points=663,full_matrix_points=816,retained_T470p=a['retained_T470p'])
    def status():
        state.update(updated=time.time(),free_bytes=shutil.disk_usage(base).free)
        h.write(root/'status.json',state)
    def step(name,cmd):
        state['phase']=name;status()
        with (root/(name+'.stdout.log')).open('x') as out,(root/(name+'.stderr.log')).open('x') as err:
            result=subprocess.run(cmd,cwd=base/'source',stdout=out,stderr=err)
        (root/(name+'.exit')).write_text(str(result.returncode)+'\n')
        if result.returncode:raise ValueError(name+' exited '+str(result.returncode))
    manifest=h.read(base/'package_hashes.json');verify_hashes(base,manifest)
    runner,importer=verify_source(base,a);verify_native(base,a)
    for name in ('audit','matrix','codespaces'):
        step('tests_'+name,[sys.executable,'-B','-X','utf8',str(base/f'source/tests/regression/test_simplemos_engineering_{name}.py')])
    step('preflight',[sys.executable,'-B','-X','utf8',str(base/'source/scripts/run_simplemos_t470p_preflight_20260929.py'),
                      '--base',str(base),'--runner',str(runner)])
    if not h.read(base/'preflight/summary.json')['passed']:raise ValueError('Preflight not qualified')
    h.write(root/'seal.json',dict(package=h.sha(base/'package_hashes.json'),assignment=h.sha(base/'assignment.json'),
        runner=h.sha(runner),importer=h.sha(importer),cases=cases,max_workers=2))
    state['phase']='running';active={}
    while state['queued'] or active:
        while state['queued'] and len(active)<2 and not state['failed']:
            if shutil.disk_usage(base).free<2*1024**3:
                state['failed'].append(dict(case=None,error='Less than 2 GiB free; dispatch stopped'));break
            try:
                verify_hashes(base,manifest)
                if h.sha(runner)!=a['runner_sha256'] or h.sha(importer)!=a['importer_sha256']:
                    raise ValueError('Frozen binary changed')
            except Exception as exc:
                state['failed'].append(dict(case=None,error=repr(exc)));break
            case=state['queued'].pop(0)
            out=(root/(case+'.stdout.log')).open('x');err=(root/(case+'.stderr.log')).open('x')
            try:
                process=subprocess.Popen([sys.executable,'-B','-X','utf8',str(base/'source/scripts/run_simplemos_engineering_matrix_20260929.py'),
                    '--base',str(base),'--case',case,'--runner',str(runner),'--importer',str(importer)],cwd=base/'source',stdout=out,stderr=err)
            except Exception as exc:
                state['failed'].append(dict(case=case,error=repr(exc)));break
            finally:out.close();err.close()
            active[case]=process;state['running'][case]=process.pid
        status()
        if not active:break
        time.sleep(10)
        for case,process in list(active.items()):
            code=process.poll()
            if code is None:continue
            (root/(case+'.exit')).write_text(str(code)+'\n')
            del active[case];del state['running'][case]
            try:
                if code:raise ValueError('Driver exited '+str(code)+'; preserved logs contain the underlying failure')
                verified_summary(base,case);state['completed'].append(case)
            except Exception as exc:state['failed'].append(dict(case=case,error=repr(exc)))
        if state['failed']:state['phase']='draining_after_failure'
    state['phase']='failed' if state['failed'] else 'completed';status()
    if state['failed']:raise ValueError('Cloud assignment stopped on preserved failure')
    collect(base,cases)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--base',required=True,type=Path);a=p.parse_args()
    try:run(a.base)
    except Exception as exc:
        h.write(a.base/'supervisor/failure.json',dict(error=repr(exc),time=time.time()))
        raise
