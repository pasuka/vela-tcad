"""Gated Release regression with integer exit codes owned by Python."""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET


def write(path,data):
    path.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')


def checked_step(name,command,cwd,out,env=None):
    with (out/(name+'.stdout.log')).open('w',encoding='utf-8') as stdout, (out/(name+'.stderr.log')).open('w',encoding='utf-8') as stderr:
        result=subprocess.run(command,cwd=cwd,env=env,stdout=stdout,stderr=stderr)
    # Preserve both successful and failed child codes; never cast a nullable
    # PowerShell Process.ExitCode to an empty string.
    (out/(name+'.exit')).write_text(str(result.returncode)+'\n')
    if result.returncode!=0:raise RuntimeError(f'{name} exited {result.returncode}')


def run(source,gate,out,toolchain):
    out.mkdir(parents=True,exist_ok=False)
    def status(phase,**extra):write(out/'status.json',dict(phase=phase,pid=os.getpid(),utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),**extra))
    try:
        status('verifying')
        for name,digest in json.loads((gate/'seal.json').read_text()).items():
            p=(gate/name).resolve()
            if not p.is_relative_to(gate.resolve()) or hashlib.sha256(p.read_bytes()).hexdigest()!=digest:raise ValueError('Gate identity mismatch')
        summary=json.loads((gate/'summary.json').read_text())
        if summary.get('numerical_passed') is not True or summary.get('points')!=816 or summary.get('curves')!=16:raise ValueError('Full numerical matrix not qualified')
        expected=json.loads((gate/'core_sha256.json').read_text())
        actual={p.relative_to(source).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for part in ('src','include') for p in (source/part).rglob('*') if p.is_file()}
        if actual!=expected:raise ValueError('Core source changed after numerical qualification')
        watched=[source/'CMakeLists.txt',source/'CMakePresets.json']
        source_seal={p.relative_to(source).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in [*watched,*[p for folder in ('tests','scripts') for p in (source/folder).rglob('*') if p.is_file() and '__pycache__' not in p.parts]]}
        write(out/'source_seal.json',dict(core=actual,regression_files=source_seal,gate_sha256=hashlib.sha256((gate/'summary.json').read_bytes()).hexdigest()))
        runner=source/'build-release/vela_example_runner.exe'
        initial=hashlib.sha256(runner.read_bytes()).hexdigest()
        env=os.environ.copy();env['PATH']=str(toolchain)+';D:/msys64/usr/bin;'+env['PATH']
        env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VELA_LINEAR_THREADS='1')
        for name,cmd in [('configure',[str(toolchain/'cmake.exe'),'--preset','windows-ucrt64-release']),('build',[str(toolchain/'cmake.exe'),'--build','--preset','windows-ucrt64-release','--parallel','2']),('inventory',[str(toolchain/'ctest.exe'),'--preset','windows-ucrt64-release','--show-only=json-v1']),('ctest',[str(toolchain/'ctest.exe'),'--preset','windows-ucrt64-release','--output-on-failure','--parallel','2','--output-junit',str(out/'ctest.xml')])]:
            status(name)
            if name=='ctest':
                inventory=json.loads((out/'inventory.stdout.log').read_text())
                if not inventory['tests']:raise ValueError('Zero selected tests')
            checked_step(name,cmd,source,out,env)
        tree=ET.parse(out/'ctest.xml').getroot()
        if int(tree.attrib.get('tests',0))<=0 or int(tree.attrib.get('failures',0)) or int(tree.attrib.get('errors',0)):raise ValueError('Invalid JUnit completion')
        if any(hashlib.sha256((source/p).read_bytes()).hexdigest()!=h for p,h in {**actual,**source_seal}.items()):raise ValueError('Source changed during regression')
        status('passed',tests=int(tree.attrib['tests']),failures=0,binary_before=initial,binary_after=hashlib.sha256(runner.read_bytes()).hexdigest())
        (out/'job.exit').write_text('0\n')
    except BaseException as exc:
        status('failed',error=repr(exc));(out/'job.exit').write_text('1\n');raise


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for k in ('source','gate','output','toolchain'):p.add_argument('--'+k,required=True,type=Path)
    a=p.parse_args();run(a.source.resolve(),a.gate.resolve(),a.output.resolve(),a.toolchain.resolve())
