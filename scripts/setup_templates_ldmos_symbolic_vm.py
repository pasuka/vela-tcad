"""Hash-checked incremental build in the existing isolated VM study checkout."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tarfile
import time


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--environment',type=Path,required=True)
    parser.add_argument('--study',type=Path,required=True)
    parser.add_argument('--resume',action='store_true')
    a=parser.parse_args();environment=a.environment.resolve();study=a.study.resolve()
    root=environment/'extrapolation_20260915';source=root/'source_a';build=root/'build-a'
    read=lambda p:json.loads(p.read_text());sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    manifest=read(root/'current_source_manifest.json');patch=read(study/'patch.json')
    # Frozen R10 binaries/results are not modified. Back up every replaced file.
    for name,digest in manifest.items():assert sha(source/name)==digest,name
    if a.resume:
        assert manifest==read(study/'source_manifest.json')
        assert not (study/'binary').exists(), 'Do not replace a frozen binary'
    else:
        (study/'prior_manifest.json').write_text(json.dumps(manifest,indent=2))
        backup=study/'source_before';backup.mkdir()
        for name in patch:
            p=source/name
            if p.exists():
                target=backup/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(str(p),str(target))
        with tarfile.open(str(study/'patch.tgz')) as archive:
            members=archive.getmembers()
            assert set(m.name for m in members)==set(patch)
            assert all(m.isfile() and not m.name.startswith('/') and '..' not in m.name.split('/') for m in members)
            archive.extractall(str(source),members=members)
        for name,digest in patch.items():
            assert sha(source/name)==digest,name
            (source/name).touch()
        manifest.update(patch)
        (root/'current_source_manifest.json').write_text(json.dumps(manifest,indent=2))
        (study/'source_manifest.json').write_text(json.dumps(manifest,indent=2))
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
    env['PATH']=str(environment/'toolchain/opt/rh/devtoolset-11/root/usr/bin')+':'+str(environment/'cmake-3.31.8-linux-x86_64/bin')+':'+env['PATH']
    env['PKG_CONFIG_PATH']=str(environment/'suitesparse_r7/pkgconfig')
    def run(argv,filename):
        path=study/filename;retry=0
        while path.exists():
            retry+=1;path=study/(filename+'.retry%d'%retry)
        with path.open('x') as log:
            start=time.perf_counter();proc=subprocess.Popen(argv,cwd=str(source),env=env,stdout=log,stderr=subprocess.STDOUT)
            while True:
                try:code=proc.wait(timeout=30);break
                except subprocess.TimeoutExpired:print(json.dumps(dict(stage=filename,elapsed=time.perf_counter()-start)),flush=True)
        assert code==0,(filename,code)
        print(json.dumps(dict(stage=filename,exit_code=code)),flush=True)
        return path
    run(['cmake','-S',str(source),'-B',str(build)],'configure.log')
    run(['cmake','--build',str(build),'--parallel','2','--target','vela_example_runner','electrothermal_probe','test_mobility','ialmob_kernel_bench'],'build.log')
    run([str(build/'test_mobility'),'[ialmob]'],'tests_mobility.log')
    runner_test=run(['ctest','--test-dir',str(build),'--output-on-failure','-R','electrothermal_dc_runner_regression'],'tests_runner.log')
    assert '1/1' in runner_test.read_text()
    flags=(build/'CMakeFiles/vela_core.dir/flags.make').read_text()
    assert '-O3' in flags and '-DNDEBUG' in flags and 'VELA_HAS_UMFPACK' in flags and '-pg' not in flags.split()
    (study/'build_flags.txt').write_text(flags)
    binary=study/'binary';binary.mkdir()
    for name in ('vela_example_runner','electrothermal_probe','ialmob_kernel_bench'):
        shutil.copy2(str(build/name),str(binary/name))
    ldd=subprocess.check_output(['ldd',str(binary/'vela_example_runner')],env=env).decode()
    assert 'not found' not in ldd
    (study/'ldd.txt').write_text(ldd)
    (study/'binary_manifest.json').write_text(json.dumps({p.name:sha(p) for p in binary.iterdir()},indent=2))
    run([str(binary/'ialmob_kernel_bench'),'2000'],'kernel_forward.json')
    run([str(binary/'ialmob_kernel_bench'),'2000','reverse'],'kernel_reverse.json')
    (study/'ready.json').write_text(json.dumps(dict(status='ready',flags=flags,binaries={p.name:sha(p) for p in binary.iterdir()}),indent=2))
    print('symbolic_release_ready',flush=True)


if __name__=='__main__':main()
